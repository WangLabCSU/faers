# Multi-group comparison and joint-group contrast.
#
# `faers_compare()` orchestrates three existing layers and adds no statistics
# of its own:
#
#   report sets      faers_reports_with()          (R/restrict-drug.R)
#   per-event counts faers_counts() / counts_db()  (R/counts.R, R/counts-db.R)
#   2x2 arithmetic   the contingency formulas of faers_phv_table()
#                                                  (R/signal.R)
#   signal metrics   phv_signal()                  (R/phv_.R)
#
# The two comparison modes differ in what the background is and must never be
# conflated, so `mode`, `comparison_id` and `background_id` travel with every
# row:
#
#   "background" - each group G against the rest of the shared scope U
#                  (a + b + c + d == |U|, one and the same U for every group)
#   "pairwise"   - two selected groups head to head
#                  (a + b + c + d == |G1| + |G2|, denominator varies per pair)
#
# Byte-parity with a hand-written single call is a hard requirement, so the
# arithmetic below reproduces R/signal.R line for line instead of re-deriving
# it, and `phv_signal()` produces every metric.

# ---------- 2x2 arithmetic (mirrors faers_phv_table) ----------

# ---------- 2x2 contingency tables ----------

# Column names used inside `data.table::[.data.table` expressions. Declared so
# that `R CMD check` does not report them as undefined globals.
utils::globalVariables(c(
    "N", "event_label", "group", "denominator", "share",
    "universe_N", "universe_share", "comparison_id", "status"
))

# Background mode: `faers_phv_table(.object = G, .full = U)`.
#   n.1 = reports of U carrying the event, n1. = |G|, n = |U|.
# The identity `G vs U` == `G vs (U \ G)` was verified against the sample data,
# which is why one cached count of U serves every group.
.compare_table_background <- function(universe_counts, group_counts, .events,
                                      n_scope, n_group) {
    data.table::setnames(universe_counts, "N", "n.1")
    data.table::setnames(group_counts, "N", "a")
    out <- merge(universe_counts, group_counts,
        by = .events, all = TRUE, allow.cartesian = TRUE
    )
    out[, a := data.table::fifelse(is.na(a), 0L, a)] # nolint
    out[, b := n_group - a] # nolint
    out[, c := n.1 - a] # nolint
    out[, d := n_scope - (n_group + n.1 - a)] # nolint
    out <- out[, !"n.1"]
    data.table::setcolorder(out, c(.events, "a", "b", "c", "d"))[]
}

# Pairwise mode: `faers_phv_table(.object = G1, .object2 = G2)`.
# Note this arithmetic carries no `a + b + c + d == const` identity of its own;
# the two groups are assumed non-overlapping, which `.overlap` guarantees.
.compare_table_pairwise <- function(counts1, counts2, .events, n1, n0) {
    data.table::setnames(counts1, "N", "a")
    data.table::setnames(counts2, "N", "c")
    out <- merge(counts1, counts2,
        by = .events, all = TRUE, allow.cartesian = TRUE
    )
    out[, c("a", "c") := lapply(.SD, function(x) { # nolint
        data.table::fifelse(is.na(x), 0L, x)
    }), .SDcols = c("a", "c")]
    out[, b := n1 - a] # nolint
    out[, d := n0 - c] # nolint
    data.table::setcolorder(out, c(.events, "a", "b", "c", "d"))[]
}

# ---------- counting ----------

# Per-group event counts, on either backend. The memory backend reuses
# `faers_keep()` + `faers_counts()`; the duckdb backend uses the TEMP-table
# join in R/compare-db.R so the user's database is never modified.
.compare_group_count <- function(.object, ids, .events, dots) {
    if (is.null(.object@db)) {
        sub <- faers_keep(.object, ids)
        return(do.call(faers_counts, c(list(.object = sub, .events = .events), dots)))
    }
    do.call(compare_counts_db, c(
        list(.object = .object, .ids = ids, .events = .events), dots
    ))
}

# Event-set restriction ------------------------------------------------
#
# `.event_set` narrows the candidate events before any counting, mirroring the
# semantics of `faers_phv_composite()` (R/signal-eventset.R): a character
# vector selects events by value, a logical vector by position, and a function
# is applied to the event column. Unlike `faers_phv_composite()` the events are
# NOT collapsed into one "Composite_Event_Set" row: the whole point of the
# event spectrum and signal matrix is one row per event, so the set only
# decides which events stay.

# Evaluate the selector against the event values and return the kept values, or
# `NULL` when no `.event_set` was given.
.compare_event_keep <- function(event_values, .event_set, .events) {
    if (is.null(.event_set)) {
        return(NULL)
    }
    if (is.function(.event_set)) {
        keep <- .event_set(event_values)
    } else if (is.logical(.event_set)) {
        if (length(.event_set) != length(event_values)) {
            cli::cli_abort(c(
                "A logical {.arg .event_set} must have one value per candidate event",
                i = "Got {length(.event_set)} value{?s} for {length(event_values)} event{?s}."
            ))
        }
        keep <- .event_set
    } else if (is.character(.event_set)) {
        keep <- event_values %chin% .event_set
    } else {
        cli::cli_abort(
            "{.arg .event_set} must be a character vector, a logical vector, or a function"
        )
    }
    if (!is.logical(keep) || length(keep) != length(event_values)) {
        cli::cli_abort(
            "{.arg .event_set} must return one logical value per candidate event"
        )
    }
    event_values[!is.na(keep) & keep]
}

# Restrict every group/scope count table to the selected events. Returning the
# empty-selection case explicitly keeps the downstream merge honest instead of
# silently comparing everything.
.compare_apply_event_set <- function(per_group, scope_counts, .event_set, .events) {
    if (is.null(.event_set)) {
        return(list(per_group = per_group, scope_counts = scope_counts, n_selected = NULL))
    }
    values <- sort(unique(c(
        scope_counts[[.events[[1L]]]],
        unlist(lapply(per_group, function(x) x[[.events[[1L]]]]), use.names = FALSE)
    )))
    keep <- .compare_event_keep(values, .event_set, .events)
    per_group <- lapply(per_group, function(x) {
        x[get(.events[[1L]]) %chin% keep]
    })
    scope_counts <- scope_counts[get(.events[[1L]]) %chin% keep]
    list(
        per_group = per_group, scope_counts = scope_counts,
        n_selected = length(keep), n_candidates = length(values)
    )
}

# Short, printable description of `.event_set` for `metadata`. The selector
# itself is not stored when it is a function or a long vector.
.compare_event_set_label <- function(.event_set) {
    if (is.null(.event_set)) {
        return("all candidate events")
    }
    if (is.function(.event_set)) {
        return("function selector")
    }
    if (is.logical(.event_set)) {
        return(sprintf("logical selector (%d TRUE)", sum(.event_set, na.rm = TRUE)))
    }
    n <- length(.event_set)
    shown <- utils::head(.event_set, 5L)
    sprintf(
        "character selector (%d value%s: %s%s)",
        n, if (n == 1L) "" else "s",
        paste(shown, collapse = ", "),
        if (n > 5L) ", ..." else ""
    )
}

# ---------- status annotation ----------

.compare_status <- function(ct) {
    out <- rep("ok", nrow(ct))
    zero <- ct$a == 0L | ct$b == 0L | ct$c == 0L | ct$d == 0L
    out[zero] <- "zero_cell"
    out[ct$a == 0L & ct$b == 0L] <- "empty_group"
    out
}

# ---------- matrices ----------

# Collapse multi-column `.events` into a single `event_label` column when the
# caller asked for more than one event column. A matrix needs one label per
# axis, so `paste(..., sep = " | ")` is used and documented; single-column
# `.events` needs no collapsing and keeps its original column name.
.compare_add_label <- function(dt, .events) {
    if (length(.events) < 2L || "event_label" %in% names(dt)) {
        return(dt)
    }
    dt[, event_label := do.call(paste, c(.SD, list(sep = " | "))), .SDcols = .events]
    dt
}

# Event x group counts into a wide matrix. `.what` selects the value column;
# for multi-column `.events` the label column is the collapsed `event_label`.
.compare_matrix <- function(profile, label_col, group_col, value_col, groups) {
    wide <- data.table::dcast(
        profile,
        stats::as.formula(sprintf("%s ~ %s", label_col, group_col)),
        value.var = value_col, fill = 0
    )
    for (g in groups) {
        if (!g %in% names(wide)) {
            data.table::set(wide, j = g, value = 0)
        }
    }
    data.table::setcolorder(wide, c(label_col, groups))
    wide[]
}

# `.top_n` truncation of the DISPLAY layer: keep the top N events per group by
# report count, take the union across groups, and collapse multi-column
# `.events` into `event_label`. Returns the truncated profile, the label column
# name to key matrices on, and the record for `metadata`.
.compare_truncate <- function(profile, .events, .top_n) {
    label_col <- if (length(.events) > 1L) "event_label" else .events[[1L]]
    profile <- .compare_add_label(profile, .events)
    if (is.null(.top_n)) {
        return(list(profile = profile, label_col = label_col, top_n = NULL, truncated = FALSE))
    }
    keep <- profile[
        order(-N, get(label_col)),
        utils::head(.SD, .top_n),
        by = "group"
    ][[label_col]]
    keep <- unique(keep)
    data.table::setattr(profile, "display_subset", TRUE)
    list(
        profile = profile[get(label_col) %chin% keep],
        label_col = label_col,
        top_n = .top_n,
        truncated = TRUE,
        n_kept = length(keep),
        rule = "per-group top-N by report count, union across groups"
    )
}

# Count / share / signal matrices. `share` and `signal` are keyed on the
# comparison actually performed, so `matrices$signal` is one column per
# `comparison_id`, while `count`/`share` are one column per group.
#
# `signal_results` is keyed on `.events`, not on the collapsed label, so the
# label is rebuilt there before the signal matrix is cast.
.compare_matrices <- function(profile, signal_results, .events, label_col,
                              .matrix_metric, group_names) {
    count <- .compare_matrix(profile, label_col, "group", "N", group_names)
    share <- .compare_matrix(profile, label_col, "group", "share", group_names)

    metric <- .matrix_metric %||% grep("_ci_low$", names(signal_results), value = TRUE)[1L]
    if (!is.null(.matrix_metric) && !.matrix_metric %in% names(signal_results)) {
        # cli interpolates from the calling frame, so bind locally.
        wanted <- .matrix_metric
        available <- setdiff(names(signal_results), .COMPARE_RESERVED_COLUMNS)
        cli::cli_abort(c(
            "{.arg .matrix_metric} {.val {wanted}} is not a column of {.field signal_results}",
            i = "Available metrics: {.val {available}}"
        ))
    }
    sig_name <- if (is.null(.matrix_metric)) "signal" else .matrix_metric
    signal_cols <- unique(signal_results$comparison_id)
    signal <- if (is.na(metric) || !length(signal_cols)) {
        data.table::data.table()
    } else {
        sub <- .compare_add_label(data.table::copy(signal_results), .events)
        sub[, comparison_id := factor(comparison_id, levels = signal_cols)]
        .compare_matrix(sub, label_col, "comparison_id", metric, signal_cols)
    }

    out <- list(count = count, share = share)
    out[[sig_name]] <- signal
    out
}

# Columns of `signal_results` that are not metrics.
.COMPARE_RESERVED_COLUMNS <- c(
    "comparison_id", "mode", "group1", "background_id", "n_total", "status"
)


#' Compare multiple drug groups against a shared background
#'
#' `faers_compare()` takes a study scope (a standardized and de-duplicated
#' [FAERSascii] object) plus one or more drug exposure definitions, builds the
#' report groups, counts every event in every group, runs disproportionality
#' analysis for each comparison and returns the results side by side: group
#' summary, event profile, the comparison definitions, the per-comparison
#' signal table, plotting-ready matrices and a metadata record.
#'
#' @details
#' ## Study scope and background
#'
#' `.object` is the scope `U`; narrow it first with [faers_restrict()]. Both
#' comparison modes only ever look inside `U`.
#'
#' * `.mode = "background"`: every group `G` is compared with the rest of `U`,
#'   i.e. the 2x2 is `a` = reports of `G` with the event, `c` = reports of
#'   `U \ G` with the event, over `n = |U|`. Each analysis uses the same `U`,
#'   reports are de-duplicated per table, and every table satisfies
#'   `a + b + c + d == |U|`. Identical to
#'   `faers_phv_table(.object = G, .full = U)`.
#' * `.mode = "pairwise"`: the background is a *selected* second group, so
#'   `n = |G1| + |G2|` and the denominator changes with every pair. Identical to
#'   `faers_phv_table(.object = G1, .object2 = G2)` once overlap is reconciled.
#'   `.comparisons` must be given explicitly; `NULL` aborts instead of silently
#'   running every pair, because per-pair denominators invite uncorrected
#'   multiple comparisons.
#'
#' The two modes are marked on every result row (`mode`, `comparison_id`,
#' `background_id`) and summarised in `comparisons`.
#'
#' ## Overlap
#'
#' `.overlap` is forwarded to [faers_group_ids()]. `"keep"` (default) is
#' correct for multi-group background comparisons: each group is compared
#' against the scope, not against another group, and the shared reports are
#' reported in `group_overlaps`. Before a *pairwise* comparison the two groups
#' must not overlap; with `"keep"` and a non-empty intersection the pair aborts
#' with the count, so pass `"exclude"` (remove the intersection from both, the
#' removed counts are recorded) or `"repartition"`/`.split = "joint"` (turn the
#' intersection into its own group).
#'
#' ## Joint-group contrast
#'
#' `faers_group_ids(.split = "joint")` splits two drugs into `A-only`,
#' `B-only`, `A+B` and `Neither` such that the four sets are mutually exclusive
#' and sum to `U`. `faers_compare()` accepts such a split either by being given
#' the four groups as `.groups` or by `.split = "joint"` directly. The result is
#' a **descriptive comparison and disproportionality analysis of co-reported
#' drugs**: it is not a formal drug-interaction test, and `Neither` exists
#' precisely so that the four-way partition stays reportable.
#'
#' ## Display subsets
#'
#' `.top_n`, when set, only truncates the *display* layer (`event_profile` and
#' `matrices`) using the rule "per group, keep the top N events by report count,
#' then take the union across groups"; ties break on the event label so the
#' result is deterministic. Truncated results are flagged with
#' `attr(., "display_subset") = TRUE` and recorded in
#' `metadata$display_subset`. `signal_results` always keeps every row so that
#' downstream filtering sees `status` instead of missing data.
#'
#' @param .object A [FAERSascii] object. Must be standardized with
#'   [faers_standardize()] and de-duplicated with [faers_dedup()].
#' @param .groups A named list of group definitions, see [faers_group_ids()].
#' @param ... Other arguments passed to [faers_counts()], e.g. `.field`,
#'   `.fn` or `.na.rm`.
#' @param .split `"none"` (default) or `"joint"`; see [faers_group_ids()].
#' @param .mode `"background"` (default) or `"pairwise"`.
#' @param .comparisons For `.mode = "pairwise"`, a list of length-2 character
#'   vectors naming the groups to contrast, e.g.
#'   `list(c("A-only", "B-only"), c("A+B", "A-only"))`. Ignored, and set to
#'   every group, in `"background"` mode.
#' @param .overlap `"keep"` (default), `"exclude"`, `"repartition"` or
#'   `"error"`; see [faers_group_ids()].
#' @param .events A character vector naming the event column(s) in `.field`,
#'   passed to [faers_counts()]. Defaults to `"pt"`. With more than one column
#'   the combination defines the event, and the matrices use a collapsed
#'   `event_label` built with `paste(..., sep = " | ")`.
#' @param .event_set `NULL` (default) to keep every candidate event, or a
#'   selector that narrows them. As in [faers_phv_composite()] it may be a
#'   character vector of event values, a logical vector with one value per
#'   candidate event, or a function applied to the event values that returns a
#'   logical vector. Unlike [faers_phv_composite()] the selected events are
#'   **not** collapsed into a single composite row: each event keeps its own row
#'   in `event_profile` and `signal_results`, which is what makes the event
#'   spectrum and the signal matrix meaningful. The candidate count, the
#'   selected count and the selector description are recorded in
#'   `metadata$event_scope`.
#' @param .min_a A whole number, `0` (default) to keep every candidate event or
#'   a reporting threshold. Comparisons whose `a` falls below it are dropped and
#'   counted in `metadata$event_scope$n_below_min_a`. Unlike
#'   [faers_phv_scan()] there is no implicit `3` here: candidate events are
#'   chosen by the caller's `.events`/`.event_set`.
#' @param .methods Methods passed to [phv_signal()]. Defaults to
#'   `c("ror", "prr")`.
#' @param .top_n `NULL` (default) for no truncation, or a whole number for the
#'   display subset described in *Details*.
#' @param .matrix_metric `NULL` to use the first available `*_ci_low` column,
#'   or the name of the metric column to place in `matrices$signal`.
#' @param .diagnose_combination Forwarded to [faers_group_ids()].
#' @param .phv_signal_params Other arguments passed to [phv_signal()].
#' @param BPPARAM Passed to [phv_signal()].
#'
#' @return A named list of [data.table][data.table::data.table] objects:
#' * `group_summary`: groups, report counts, share of scope, overlap removals,
#'   definitions and the single-row combination count.
#' * `event_profile`: per group and event, `N` reports, the `denominator` and
#'   the within-group `share`, plus the same quantities over the whole scope
#'   (`universe_N`, `universe_share`). Every group is laid out over the scope's
#'   candidate events, so absent combinations appear as `N = 0` and the table is
#'   a complete events x groups grid that lines up with `signal_results`.
#' * `comparisons`: which group was compared with which background, the group
#'   sizes and the overlap removal counts.
#' * `signal_results`: one row per comparison x event with `a`, `b`, `c`, `d`,
#'   the `expected` count, every metric from [phv_signal()] and a `status`
#'   column (`"ok"`, `"zero_cell"`, `"empty_group"`).
#' * `matrices`: `count`, `share` and one signal matrix (named after
#'   `.matrix_metric`, or `signal` for the default `*_ci_low` column), each
#'   event x group/comparison.
#' * `metadata`: scope, groups, events, comparison and denominator definitions,
#'   versions, missing-data and truncation records.
#' * `group_overlaps`: pairwise intersections of the exposure definitions.
#'
#' @examples
#' data <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
#' data <- faers_dedup(data)
#'
#' # three drug groups, each against the same scope
#' res <- suppressMessages(faers_compare(
#'     data,
#'     list(
#'         insulin = "insulin|humalog",
#'         biguanide = "metformin",
#'         sulfonylurea = "glyburide|glipizide"
#'     ),
#'     .events = "pt", .min_a = 2L
#' ))
#' res$group_summary
#' res$comparisons
#' head(res$signal_results)
#' res$matrices$count[1:5, 1:4]
#'
#' # two drugs crossed into four mutually exclusive groups
#' joint <- suppressMessages(faers_compare(
#'     data, list(A = "insulin|humalog", B = "metformin"),
#'     .split = "joint", .mode = "pairwise",
#'     .comparisons = list(c("A-only", "B-only")),
#'     .events = "pt", .min_a = 1L
#' ))
#' joint$comparisons
#'
#' @seealso [faers_group_ids()], [faers_phv_table()], [faers_phv_signal()],
#'   [phv_signal()], [faers_phv_scan()], [faers_reports_with()],
#'   [faers_restrict()].
#' @export
faers_compare <- function(.object,
                          .groups,
                          ...,
                          .split = c("none", "joint"),
                          .mode = c("background", "pairwise"),
                          .comparisons = NULL,
                          .overlap = c("keep", "exclude", "repartition", "error"),
                          .events = "pt",
                          .event_set = NULL,
                          .min_a = 0L,
                          .methods = c("ror", "prr"),
                          .top_n = NULL,
                          .matrix_metric = NULL,
                          .diagnose_combination = TRUE,
                          .phv_signal_params = list(),
                          BPPARAM = BiocParallel::SerialParam()) {
    .split <- match.arg(.split)
    .mode <- match.arg(.mode)
    .overlap <- match.arg(.overlap)
    assert_s3_class(.object, "FAERSascii")
    assert_character(.events, allow_na = FALSE)
    if (!length(.events) || any(!nzchar(.events))) {
        cli::cli_abort("{.arg .events} must be a non-empty character vector")
    }
    assert_number_whole(.min_a, min = 0)
    assert_number_whole(.top_n, min = 1, allow_null = TRUE)
    assert_string(.matrix_metric, allow_null = TRUE, allow_empty = FALSE)
    assert_(.phv_signal_params, is.list, "a list")
    assert_inclusive(.methods, SCAN_ALL_METHODS, allow_null = TRUE)
    if (!.object@standardization) {
        cli::cli_abort("{.arg .object} must be standardized using {.fn faers_standardize}")
    }
    if (!.object@deduplication) {
        cli::cli_abort(c(
            "{.arg .object} must be de-duplicated using {.fn faers_dedup}",
            i = "Comparing un-deduplicated groups would bias every contingency table"
        ))
    }
    .events_check <- .compare_check_events(.object, .events)

    # ---- (1) report groups ----
    groups <- faers_group_ids(
        .object, .groups,
        .split = .split, .overlap = .overlap,
        .diagnose_combination = .diagnose_combination
    )
    ids <- groups$ids
    scope_ids <- groups$scope_ids
    n_scope <- length(scope_ids)

    # ---- (2) the comparison plan ----
    plan <- .compare_plan(.mode, .comparisons, names(ids))

    # ---- (3) counts (scope counts are computed once and reused) ----
    dots <- list(...)
    scope_counts <- .compare_group_count(.object, scope_ids, .events, dots)
    per_group <- lapply(ids, function(x) {
        .compare_group_count(.object, x, .events, dots)
    })

    # ---- (3b) optional candidate-event restriction ----
    sel <- .compare_apply_event_set(per_group, scope_counts, .event_set, .events)
    per_group <- sel$per_group
    scope_counts <- sel$scope_counts

    # ---- (4) contingency tables per comparison ----
    tables <- .compare_tables(
        plan, per_group, scope_counts, ids, .events, n_scope, .object
    )

    # ---- (5) statistics ----
    signal <- .compare_signal(tables, .methods, .phv_signal_params, BPPARAM)

    # ---- (6) results ----
    group_summary <- groups$summary
    signal_results <- .compare_signal_results(signal, plan, ids, .events, .min_a)
    comparisons <- .compare_comparison_table(
        plan, ids, n_scope, .overlap, groups$raw_ids
    )
    event_profile <- .compare_event_profile(
        groups, per_group, scope_counts, .events, n_scope,
        .field = dots$.field %||% "reac"
    )
    # `.top_n` truncates the display layer only; `signal_results` keeps every row.
    trunc <- .compare_truncate(event_profile, .events, .top_n)
    matrices <- .compare_matrices(
        trunc$profile, signal_results, .events, trunc$label_col,
        .matrix_metric, names(ids)
    )
    metadata <- .compare_metadata(
        .object, .groups, groups, plan, .mode, .split, .overlap,
        .events, .event_set, .min_a, .methods, .phv_signal_params,
        .top_n, .matrix_metric,
        dots$.field %||% "reac", n_scope, scope_counts, signal_results,
        .events_check, trunc, sel
    )

    list(
        group_summary = group_summary,
        event_profile = trunc$profile,
        comparisons = comparisons,
        signal_results = signal_results,
        matrices = matrices,
        metadata = metadata,
        group_overlaps = groups$overlaps
    )
}

# ---------- validation ----------

# Event columns must exist in the standardized `reac` data, as in
# `faers_phv_scan()`. Skipped for a non-`reac` `.field`, where the caller is on
# their own (matching `faers_counts()`'s behaviour).
.compare_check_events <- function(.object, .events, .field = "reac") {
    if (!identical(.field, "reac")) {
        return(list(field = .field, allowed = NA_character_))
    }
    allowed <- c(
        "pt", "meddra_code", "meddra_pt",
        names(faers_meddra(.object, use = "hierarchy"))
    )
    bad <- setdiff(.events, allowed)
    if (length(bad)) {
        cli::cli_abort(c(
            "{.arg .events} must be column{?s} of the standardized {.field reac} data",
            x = "Unknown: {.val {bad}}",
            i = "Available columns: {.val {allowed}}"
        ))
    }
    list(field = .field, allowed = allowed)
}

# ---------- comparison plan ----------

# What is compared against what. The two modes are kept apart here so no
# downstream step has to guess: a comparison is always an explicit list.
.compare_plan <- function(.mode, .comparisons, group_names) {
    if (identical(.mode, "background")) {
        if (!is.null(.comparisons)) {
            cli::cli_inform(c(
                "{.arg .comparisons} is ignored with {.val .mode = \"background\"}",
                i = "Every group is compared against the rest of the study scope"
            ))
        }
        return(lapply(group_names, function(g) {
            list(
                comparison_id = paste0(g, "_vs_rest_of_scope"),
                mode = "background",
                group1 = g, group2 = NA_character_,
                background_id = "rest_of_scope"
            )
        }))
    }
    # ---- pairwise ----
    if (is.null(.comparisons)) {
        cli::cli_abort(c(
            "{.arg .comparisons} must be supplied with {.val .mode = \"pairwise\"}",
            i = "Each pair is contrasted against a different background, so they are never run implicitly",
            i = "Available groups: {.val {group_names}}",
            i = "Example: {.code .comparisons = list(c(\"{group_names[[1L]]}\", \"{group_names[[length(group_names)]]}\"))}"
        ))
    }
    if (!is.list(.comparisons) || !length(.comparisons)) {
        cli::cli_abort("{.arg .comparisons} must be a non-empty list of length-2 character vectors")
    }
    lapply(seq_along(.comparisons), function(i) {
        pair <- .comparisons[[i]]
        if (!is.character(pair) || length(pair) != 2L || anyNA(pair)) {
            cli::cli_abort(c(
                "{.arg .comparisons}[[{i}]] must be a character vector of length 2",
                i = "Available groups: {.val {group_names}}"
            ))
        }
        bad <- setdiff(pair, group_names)
        if (length(bad)) {
            cli::cli_abort(c(
                "{.arg .comparisons}[[{i}]] refers to unknown group{?s} {.val {bad}}",
                i = "Available groups: {.val {group_names}}"
            ))
        }
        if (identical(pair[[1L]], pair[[2L]])) {
            cli::cli_abort("{.arg .comparisons}[[{i}]] names the same group twice")
        }
        list(
            comparison_id = paste0(pair[[1L]], "_vs_", pair[[2L]]),
            mode = "pairwise",
            group1 = pair[[1L]], group2 = pair[[2L]],
            background_id = pair[[2L]]
        )
    })
}

# ---------- contingency tables ----------

# One 2x2 per comparison. Returns a list aligned with `plan`, each element
# carrying its plan entry, the table, and the `n_total` the table conserves.
.compare_tables <- function(plan, per_group, scope_counts, ids, .events,
                            n_scope, .object) {
    lapply(plan, function(cmp) {
        if (identical(cmp$mode, "background")) {
            ct <- .compare_table_background(
                data.table::copy(scope_counts),
                data.table::copy(per_group[[cmp$group1]]),
                .events, n_scope, length(ids[[cmp$group1]])
            )
            return(list(plan = cmp, ct = ct, n_total = n_scope))
        }
        # Pairwise: the arithmetic assumes disjoint groups, so refuse to build a
        # table whose denominator would silently double-count shared reports.
        overlap <- intersect(ids[[cmp$group1]], ids[[cmp$group2]])
        if (length(overlap)) {
            cli::cli_abort(c(
                "Groups {.val {cmp$group1}} and {.val {cmp$group2}} share {length(overlap)} report{?s}",
                i = "A pairwise table needs disjoint groups: use {.val .overlap = \"exclude\"}, {.val \"repartition\"} or {.code .split = \"joint\"}"
            ))
        }
        n1 <- length(ids[[cmp$group1]])
        n0 <- length(ids[[cmp$group2]])
        ct <- .compare_table_pairwise(
            data.table::copy(per_group[[cmp$group1]]),
            data.table::copy(per_group[[cmp$group2]]),
            .events, n1, n0
        )
        list(plan = cmp, ct = ct, n_total = n1 + n0)
    })
}

# ---------- statistics ----------

# Run `phv_signal()` once per comparison. Statistics are never reimplemented:
# the table is handed to the same function a single-group analysis would use.
# This holds for an empty table too - `phv_signal()` returns its full schema on
# zero-length input - so an empty comparison produces the same columns as a
# populated one and no result shape depends on the data.
.compare_signal <- function(tables, .methods, .phv_signal_params, BPPARAM) {
    lapply(tables, function(tbl) {
        ct <- tbl$ct
        sig <- do.call(
            phv_signal,
            c(
                ct[, c("a", "b", "c", "d")],
                list(methods = .methods, BPPARAM = BPPARAM),
                .phv_signal_params
            )
        )
        list(plan = tbl$plan, ct = ct, n_total = tbl$n_total, signal = sig)
    })
}

# ---------- result tables ----------

# One row per comparison x event: the contingency table, the metrics, the
# provenance (`mode`, `group1`, `background_id`, `n_total`) and the `status`
# flag. `.min_a` drops rows below the reporting threshold and reports the count
# through `attr(, "n_below_min_a")` rather than hiding it.
.compare_signal_results <- function(signal, plan, ids, .events, .min_a) {
    n_below <- 0L
    out <- data.table::rbindlist(lapply(signal, function(tbl) {
        ct <- tbl$ct
        if (!nrow(ct)) {
            return(NULL)
        }
        cmp <- tbl$plan
        sig <- tbl$signal
        res <- data.table::data.table(
            comparison_id = cmp$comparison_id,
            mode = cmp$mode,
            group1 = cmp$group1,
            background_id = cmp$background_id,
            n_total = as.integer(tbl$n_total)
        )
        res <- cbind(res, ct, sig)
        res
    }), use.names = TRUE, fill = TRUE)

    if (!nrow(out)) {
        # Every comparison had an empty contingency table. Rebuild the schema
        # from the (zero-row) signal tables so an empty result carries exactly
        # the columns a populated one does, metric columns included - the
        # metric names come from `phv_signal()` instead of being hardcoded, so
        # e.g. `chisq` is never returned without its `chisq_pvalue` companion.
        sig_names <- if (length(signal)) names(signal[[1L]]$signal) else character()
        out <- data.table::data.table(
            comparison_id = character(), mode = character(), group1 = character(),
            background_id = character(), n_total = integer()
        )
        for (e in .events) data.table::set(out, j = e, value = character())
        for (v in c("a", "b", "c", "d")) data.table::set(out, j = v, value = integer())
        for (v in sig_names) data.table::set(out, j = v, value = numeric())
        data.table::set(out, j = "status", value = character())
        data.table::setattr(out, "n_below_min_a", 0L)
        return(out[])
    }

    below <- out$a < .min_a
    n_below <- sum(below)
    out[, status := .compare_status(.SD), .SDcols = c("a", "b", "c", "d")]
    if (n_below) {
        out <- out[!below]
    }
    data.table::setattr(out, "n_below_min_a", as.integer(n_below))
    out[]
}

# Which group was compared against which background, with the sizes and the
# overlap removals that produced them.
#
# `n_overlap_removed` is the number of reports the overlap strategy took out of
# the two compared groups before the table was built. It is measured on the RAW
# exposure sets, because by the time the reconciled `ids` exist the shared
# reports are already gone from both groups (`"exclude"`) or promoted into their
# own group (`"repartition"`). With `"keep"` nothing is removed, and a shared
# report would be counted in both arms -- which is exactly what `"keep"` means.
.compare_comparison_table <- function(plan, ids, n_scope, .overlap,
                                      raw_ids = NULL, raw_overlap = NULL) {
    data.table::rbindlist(lapply(plan, function(cmp) {
        if (identical(cmp$mode, "background")) {
            return(data.table::data.table(
                comparison_id = cmp$comparison_id,
                mode = cmp$mode,
                group1 = cmp$group1,
                group2 = NA_character_,
                background_id = cmp$background_id,
                n1 = as.integer(length(ids[[cmp$group1]])),
                n0 = as.integer(n_scope - length(ids[[cmp$group1]])),
                n_total = as.integer(n_scope),
                n_overlap_removed = 0L,
                n_overlap_shared = 0L,
                strategy = .overlap
            ))
        }
        raw_shared <- 0L
        if (!is.null(raw_ids)) {
            raw_shared <- length(intersect(raw_ids[[cmp$group1]], raw_ids[[cmp$group2]]))
        }
        removed <- if (.overlap %chin% c("exclude", "repartition")) raw_shared else 0L
        data.table::data.table(
            comparison_id = cmp$comparison_id,
            mode = cmp$mode,
            group1 = cmp$group1,
            group2 = cmp$group2,
            background_id = cmp$background_id,
            n1 = as.integer(length(ids[[cmp$group1]])),
            n0 = as.integer(length(ids[[cmp$group2]])),
            n_total = as.integer(length(ids[[cmp$group1]]) + length(ids[[cmp$group2]])),
            n_overlap_removed = as.integer(removed),
            n_overlap_shared = as.integer(raw_shared),
            strategy = .overlap
        )
    }), fill = TRUE)
}

# Per group x event: the group's report count, its denominator (the group's own
# report count, so `share` is a within-group proportion) and the same quantities
# over the whole scope.
#
# Every group is laid out over the SCOPE's candidate events, with `N = 0` where
# the group has no report for an event. That is what makes `event_profile` and
# the count/share matrices a complete events x groups grid consistent with
# `signal_results`, whose rows also come from the scope counts.
.compare_event_profile <- function(groups, per_group, scope_counts, .events,
                                   n_scope, .field = "reac") {
    scope <- data.table::copy(scope_counts)
    data.table::setnames(scope, "N", "universe_N")
    cn <- .events
    out <- data.table::rbindlist(lapply(names(per_group), function(nm) {
        g <- data.table::copy(per_group[[nm]])
        grid <- merge(
            scope[, c(cn, "universe_N"), with = FALSE], g,
            by = cn, all.x = TRUE, sort = FALSE
        )
        grid[is.na(N), N := 0L]
        grid[, group := nm]
        grid
    }), use.names = TRUE, fill = TRUE)
    if (!nrow(out)) {
        # Keep the full schema so an empty `.event_set` selection still yields a
        # well-formed, zero-row profile instead of a table whose columns are
        # missing for every downstream consumer.
        empty <- data.table::data.table(group = character())
        for (e in .events) data.table::set(empty, j = e, value = character())
        for (v in c("N", "denominator")) data.table::set(empty, j = v, value = integer())
        for (v in c("share", "universe_N", "universe_share")) {
            data.table::set(empty, j = v, value = numeric())
        }
        data.table::setcolorder(
            empty,
            c("group", .events, "N", "denominator", "share", "universe_N", "universe_share")
        )
        return(empty[])
    }
    out[, denominator := as.integer(length(groups$ids[[.BY$group]])), by = "group"]
    out[, share := N / denominator]
    out[, universe_share := universe_N / n_scope]
    data.table::setcolorder(
        out,
        c("group", .events, "N", "denominator", "share", "universe_N", "universe_share")
    )
    data.table::setorderv(
        out,
        c("group", "N", .events),
        order = c(1L, -1L, rep(1L, length(.events)))
    )
    out[]
}

# ---------- metadata ----------

# Everything a reader needs to reproduce the run: what the scope was, how the
# groups were defined, which comparison each table belongs to, which rules
# filtered rows, and which software versions produced it.
.compare_metadata <- function(.object, raw_groups, groups, plan, .mode, .split,
                              .overlap, .events, .event_set, .min_a, .methods,
                              .phv_signal_params, .top_n, .matrix_metric,
                              .field, n_scope, scope_counts, signal_results,
                              .events_check, trunc, sel) {
    meddra <- faers_meddra(.object, use = "hierarchy")
    list(
        call = match.call(expand.dots = FALSE),
        faers_version = as.character(utils::packageVersion("faers")),
        r_version = paste(R.version$major, R.version$minor, sep = "."),
        timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
        mode = .mode,
        split = .split,
        overlap_strategy = .overlap,
        groups = raw_groups,
        scope = list(
            n_reports_universe = as.integer(n_scope),
            year = .object@year,
            quarter = .object@quarter,
            standardization = .object@standardization,
            deduplication = .object@deduplication,
            n_deleted_cases = length(faers_deleted_cases(.object)),
            meddra_version = if (is.null(meddra)) NA_character_ else meddra_version(.object@meddra),
            n_meddra_hierarchy_columns = if (is.null(meddra)) 0L else ncol(meddra)
        ),
        events = .events,
        event_scope = list(
            field = .field,
            allowed_columns = .events_check$allowed,
            event_set = .compare_event_set_label(.event_set),
            n_candidate_events = sel$n_candidates %||% nrow(scope_counts),
            n_selected_events = sel$n_selected %||% nrow(scope_counts),
            n_universe_events = nrow(scope_counts),
            min_a = .min_a,
            n_below_min_a = attr(signal_results, "n_below_min_a") %||% 0L
        ),
        methods = .methods %||% SCAN_ALL_METHODS,
        phv_signal_params = .phv_signal_params,
        denominator_definition = switch(.mode,
            background = "n = |U|: each group is compared against the rest of the study scope, a + b + c + d == |U|",
            pairwise = "n = |G1| + |G2|: the background is the selected second group, a + b + c + d == |G1| + |G2|"
        ),
        background_definition = switch(.mode,
            background = "rest_of_scope = U minus the compared group; the same U is used for every group",
            pairwise = "the other group of the pair; the denominator changes with every comparison"
        ),
        comparisons = vapply(plan, `[[`, character(1L), "comparison_id"),
        overlap_handling = list(
            strategy = .overlap,
            n_shared_pairs = if (is.null(groups$overlaps)) 0L else nrow(groups$overlaps),
            n_reports_shared = if (is.null(groups$overlaps)) 0L else
                as.integer(sum(groups$overlaps$n_intersection)),
            pairwise = .overlap %chin% c("exclude", "repartition", "error"),
            note = switch(.overlap,
                keep = paste(
                    "Overlapping reports are kept in every group. In background",
                    "mode each group still uses the same scope, so no report is",
                    "counted twice in one table; in pairwise mode a shared report",
                    "appears in both arms."
                ),
                exclude = paste(
                    "Reports shared by two compared groups are removed from both",
                    "before the table is built; the removed counts are reported in",
                    "`comparisons$n_overlap_removed`."
                ),
                repartition = paste(
                    "Reports shared by two compared groups are promoted into an",
                    "extra, mutually exclusive group (`<A>&<B>`) instead of being",
                    "dropped."
                ),
                error = "An overlapping pair aborts the run."
            )
        ),
        missing = list(
            na_rm = FALSE,
            note = "NA event values are kept as their own row, matching faers_counts() defaults"
        ),
        truncation = list(
            top_n = .top_n,
            matrices_applied = trunc$truncated,
            event_profile_applied = trunc$truncated,
            signal_results_applied = FALSE,
            rule = trunc$rule %||% NA_character_,
            n_events_kept = trunc$n_kept %||% NA_integer_
        ),
        display_subset = list(
            event_profile = trunc$truncated,
            matrices = trunc$truncated,
            signal_results = FALSE,
            top_n = .top_n,
            rule = "per-group top-N by report count, union across groups"
        ),
        matrix_metric = .matrix_metric %||% "first *_ci_low column",
        interaction_test = "None. Joint groups are a descriptive comparison and disproportionality analysis of co-reported drugs; a formal interaction test requires its own method definition."
    )
}

