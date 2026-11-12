# Report-group construction for multi-group comparison.
#
# A *report group* is a named exposure definition: a set of drug name patterns
# (regular expressions, case-insensitive), optionally a per-drug `role_cod`
# restriction, an AND/OR combination rule and a set of report-level exclusion
# patterns. The definitions are resolved with `faers_reports_with()` (see
# R/restrict-drug.R), so the matching semantics - including the rule that under
# `.logic = "and"` two patterns MAY be satisfied by the SAME drug row
# (fixed-dose combinations / multi-ingredient `prod_ai` entries) and that
# `.exclude_drugs` is report-level and role-agnostic - are inherited verbatim
# rather than reimplemented.
#
# This file only builds and reconciles *report sets*:
#   `faers_group_ids()` -> ids + membership + summary + pairwise overlaps
#
# The contingency tables, statistics and matrices live in R/compare.R.

# ---------- group definition validation ----------

# Columns referenced inside `data.table::[.data.table` expressions.
utils::globalVariables(c("n_intersection"))

.COMPARE_GROUP_KEYS <- c(".drugs", ".roles", ".logic", ".exclude_drugs", ".drug_col")
.COMPARE_OVERLAP_STRATEGIES <- c("keep", "exclude", "repartition", "error")

# Normalize `.groups` into a named list of fully-specified definitions. A group
# may be given as a plain character vector (shorthand for "these patterns, no
# role restriction, AND, no exclusion") or as a list of the keys above.
.compare_check_groups <- function(.groups, .split) {
    if (!is.list(.groups) || !length(.groups)) {
        cli::cli_abort("{.arg .groups} must be a non-empty named list of group definitions")
    }
    nms <- names(.groups)
    if (is.null(nms) || anyNA(nms) || any(!nzchar(nms))) {
        cli::cli_abort("{.arg .groups} must be fully named, e.g. {.code list(A = \"warfarin\", B = \"metformin\")}")
    }
    if ("primaryid" %in% nms) {
        cli::cli_abort("{.arg .groups} cannot use the reserved name {.val primaryid}; it identifies reports in {.field membership}")
    }
    dup <- unique(nms[duplicated(nms)])
    if (length(dup)) {
        cli::cli_abort(c(
            "{.arg .groups} contains duplicated group names",
            x = "Duplicated: {.val {dup}}"
        ))
    }
    # `cli` interpolates `{.val {x}}` from the calling frame, which does not see
    # namespace constants, so the allowed keys are bound locally first.
    allowed_keys <- .COMPARE_GROUP_KEYS
    out <- lapply(seq_along(.groups), function(i) {
        nm <- nms[[i]]
        def <- .groups[[i]]
        if (is.character(def)) {
            def <- list(.drugs = def)
        }
        if (!is.list(def)) {
            cli::cli_abort(c(
                "Group {.val {nm}} must be a character vector of drug patterns or a named list",
                i = "Allowed keys: {.val {allowed_keys}}"
            ))
        }
        unknown <- setdiff(names(def), .COMPARE_GROUP_KEYS)
        if (length(unknown)) {
            cli::cli_abort(c(
                "Group {.val {nm}} uses unknown key{?s} {.val {unknown}}",
                i = "Allowed keys: {.val {allowed_keys}}"
            ))
        }
        drugs <- .reports_with_check_drugs(def$.drugs)
        roles <- .reports_with_check_roles(def$.roles, drugs)
        logic <- def$.logic %||% "and"
        logic <- restrict_check_values(logic, c("and", "or"), sprintf('.groups[["%s"]][[".logic"]]', nm))
        if (length(logic) != 1L) {
            cli::cli_abort("Group {.val {nm}}: {.arg .logic} must be a single string, {.val \"and\"} or {.val \"or\"}")
        }
        list(
            name = nm,
            .drugs = drugs,
            .roles = roles,
            .logic = logic,
            .exclude_drugs = .reports_with_check_exclude(def$.exclude_drugs),
            .drug_col = def$.drug_col
        )
    })
    names(out) <- nms
    if (identical(.split, "joint") && length(out) != 2L) {
        cli::cli_abort(c(
            "{.arg .split = \"joint\"} requires exactly 2 groups, got {length(out)}",
            i = "Supply the two drugs to be crossed, e.g. {.code list(A = \"warfarin\", B = \"metformin\")}"
        ))
    }
    out
}

# Resolve one group definition into a sorted unique `primaryid` vector.
# `faers_reports_with()` emits its own per-pattern message; silence it here so a
# multi-group run does not print one message per group and per comparison.
.compare_group_ids <- function(.object, def) {
    # `faers_reports_with()` accepts only `NULL` for "no role restriction", not
    # an empty list, so normalize back before delegating.
    roles <- if (length(def$.roles)) def$.roles else NULL
    excl <- if (length(def$.exclude_drugs)) def$.exclude_drugs else NULL
    suppressMessages(faers_reports_with(
        .object,
        .drugs = def$.drugs,
        .roles = roles,
        .logic = def$.logic,
        .exclude_drugs = excl,
        .drug_col = def$.drug_col,
        .return = "ids"
    ))
}

# The drug column a definition actually matched against, after the
# `drugname` -> `prod_ai` fallback that `faers_reports_with()` performs.
.compare_group_drug_col <- function(.object, def) {
    .reports_with_pick_col(.object, def$.drug_col)
}

# ---------- joint split ----------

.COMPARE_ONLY <- "-only"
.COMPARE_BOTH <- "+"
.COMPARE_NEITHER <- "Neither"

# Split two exposure sets into the four mutually exclusive joint groups. The
# complement is taken against the SAME exposure definition (`B_ids`), never
# against "any drug row mentioning B": with a role-restricted B only its
# role-qualified rows build `B_ids`, so a report whose B record carries another
# role is correctly absent from B and its complement.
.compare_split_joint <- function(ids, scope_ids) {
    a <- ids[[1L]]
    b <- ids[[2L]]
    nms <- names(ids)
    out <- list(
        setdiff(a, b),
        setdiff(b, a),
        intersect(a, b),
        setdiff(scope_ids, union(a, b))
    )
    names(out) <- c(
        paste0(nms[[1L]], .COMPARE_ONLY),
        paste0(nms[[2L]], .COMPARE_ONLY),
        paste0(nms[[1L]], .COMPARE_BOTH, nms[[2L]]),
        .COMPARE_NEITHER
    )
    out
}

# ---------- overlap reconciliation ----------

# All pairwise intersections, for diagnosis and for the comparison layer.
.compare_overlap_table <- function(ids, strategy) {
    nms <- names(ids)
    if (length(nms) < 2L) {
        return(data.table::data.table(
            group1 = character(), group2 = character(),
            n_intersection = integer(), strategy = character()
        ))
    }
    pairs <- utils::combn(nms, 2L, simplify = FALSE)
    out <- data.table::rbindlist(lapply(pairs, function(p) {
        n <- length(intersect(ids[[p[[1L]]]], ids[[p[[2L]]]]))
        data.table::data.table(
            group1 = p[[1L]], group2 = p[[2L]],
            n_intersection = as.integer(n), strategy = strategy
        )
    }))
    out[]
}

# Globally mutually exclusive sets: each group loses every report listed in any
# other group. Used to turn overlapping exposure definitions into a partition.
.compare_overlap_exclude <- function(ids) {
    nms <- names(ids)
    out <- lapply(nms, function(nm) {
        others <- setdiff(nms, nm)
        setdiff(ids[[nm]], Reduce(union, ids[others], init = character()))
    })
    names(out) <- nms
    out
}

# Turn the intersection of exactly two overlapping groups into a third group so
# nothing is dropped: `<A>&<B>` joins `A-only` and `B-only`. Restricted to two
# groups because the partition of three or more overlapping sets is ambiguous;
# use `.split = "joint"` for a two-drug cross instead.
.compare_overlap_repartition <- function(ids) {
    nms <- names(ids)
    if (length(nms) != 2L) {
        cli::cli_abort(c(
            "{.val .overlap = \"repartition\"} is only defined for exactly 2 groups",
            i = "For a two-drug cross use {.code .split = \"joint\"}; otherwise reconcile overlaps pairwise"
        ))
    }
    a <- ids[[1L]]
    b <- ids[[2L]]
    both <- intersect(a, b)
    out <- list(setdiff(a, b), setdiff(b, a), both)
    names(out) <- c(
        paste0(nms[[1L]], .COMPARE_ONLY),
        paste0(nms[[2L]], .COMPARE_ONLY),
        paste0(nms[[1L]], "&", nms[[2L]])
    )
    out
}

.compare_apply_overlap <- function(ids, strategy) {
    switch(strategy,
        keep = ids,
        exclude = .compare_overlap_exclude(ids),
        repartition = .compare_overlap_repartition(ids),
        error = {
            ov <- .compare_overlap_table(ids, strategy)
            clash <- ov[n_intersection > 0L]
            if (nrow(clash)) {
                cli::cli_abort(c(
                    "Overlapping report groups with {.val .overlap = \"error\"}",
                    x = "{clash$group1} / {clash$group2}: {clash$n_intersection} shared report{?s}",
                    i = "Use {.val .overlap = \"exclude\"} or {.val \"repartition\"} to reconcile"
                ))
            }
            ids
        },
        cli::cli_abort("Unknown overlap strategy {.val {strategy}}")
    )
}

# ---------- single-row combination diagnosis ----------

# Row-level group coverage: for every drug row, does it alone satisfy the whole
# definition? Under `.logic = "and"` that means every pattern (respecting each
# pattern's role restriction), under `"or"` at least one.
.compare_row_covers <- function(drug, def) {
    if (!nrow(drug)) {
        return(logical())
    }
    hits <- lapply(def$.drugs, function(p) {
        .reports_with_row_hit(drug, p, def$.roles[[p]], def$.drug_col)
    })
    hits <- do.call(cbind, hits)
    if (identical(def$.logic, "and")) {
        rowSums(hits) == ncol(hits)
    } else {
        rowSums(hits) >= 1L
    }
}

# How many reports in `ids_both` have ONE drug row that satisfies both
# definitions at once - i.e. how much of the joint group comes from a
# fixed-dose combination product rather than from two separately recorded
# exposures. This is a description of the exposure data, not a statistical
# interaction test, and it is only computed when `.diagnose_combination`.
.compare_single_row_both <- function(.object, defs, ids_both) {
    if (!length(ids_both)) {
        return(0L)
    }
    resolved <- lapply(defs, function(d) {
        d$.drug_col <- .compare_group_drug_col(.object, d)
        d
    })
    cols <- unique(c(vapply(resolved, `[[`, character(1L), ".drug_col"), "role_cod"))
    drug <- compare_drug_rows(.object, cols)
    keep_ids <- unique(as.character(ids_both))
    drug <- drug[primaryid %chin% keep_ids]
    if (!nrow(drug)) {
        return(0L)
    }
    covers <- lapply(resolved, function(d) .compare_row_covers(drug, d))
    both <- Reduce(`&`, covers)
    as.integer(length(unique(drug$primaryid[both])))
}

# ---------- main entry ----------

#' Build report groups from drug exposure definitions
#'
#' `faers_group_ids()` resolves one or more named drug exposure definitions into
#' report (`primaryid`) sets within a fixed study scope, diagnoses how much the
#' groups overlap, and optionally splits two drugs into the four mutually
#' exclusive joint groups `A-only`, `B-only`, `A+B` and `Neither`.
#'
#' @details
#' ## Where the study scope comes from
#'
#' The scope `U` is `.object` itself: narrow it first with [faers_restrict()]
#' (period, demographic or geographic conditions) so that every group is built
#' inside one and the same report universe. No second "full database" argument
#' exists - [faers_compare()] compares each group against the rest of `U`.
#'
#' ## Group definitions
#'
#' Each element of `.groups` is either a character vector of drug name patterns
#' or a list with the keys `.drugs`, `.roles`, `.logic`, `.exclude_drugs` and
#' `.drug_col`, all of which carry the meaning documented in
#' [faers_reports_with()] because that function does the matching. In
#' particular matching is a case-insensitive regular expression, `.roles`
#' restricts which `role_cod` values may satisfy a pattern, and under
#' `.logic = "and"` two patterns **may be satisfied by the same drug row**
#' (fixed-dose combination products). A report matching both an `A` and a `B`
#' definition is therefore placed in the joint group `A+B`; the
#' `n_single_row_both` column of `group_summary` counts how many joint-group
#' reports come from one such multi-ingredient row.
#'
#' ## Joint split
#'
#' With `.split = "joint"` (exactly two groups) the four sets are
#' `A-only = A \ B`, `B-only = B \ A`, `A+B = A n B` and
#' `Neither = U \ (A u B)`. The complement is taken against the **same**
#' exposure definition that built the other group, so a role-restricted `B`
#' (`PS`/`SS` only) contributes only its role-qualified reports to `B` - a
#' report whose B record carries another role is not silently treated as
#' B-exposed. The four sets are mutually exclusive and sum to `U`.
#'
#' ## `.overlap`
#'
#' `"keep"` (default) returns the exposure sets as defined and only reports the
#' pairwise intersections in `group_overlaps`; this is what multi-drug
#' background comparisons need, since each group is compared against the shared
#' scope rather than against another group. `"exclude"` makes the sets mutually
#' exclusive by removing, from every group, all reports listed in any other
#' group. `"repartition"` (two groups only) keeps the information by turning the
#' intersection into a third group `<A>&<B>` alongside `A-only` and `B-only`.
#' `"error"` aborts when any two groups share a report.
#'
#' @param .object A [FAERSascii] object, standardized with
#'   [faers_standardize()]. It defines the study scope `U`.
#' @param .groups A named list of group definitions, see *Details*.
#'   The name `primaryid` is reserved for the report ID column.
#' @param .split `"none"` (default) keeps the groups as given; `"joint"` splits
#'   exactly two groups into `A-only`, `B-only`, `A+B` and `Neither`.
#' @param .overlap `"keep"`, `"exclude"`, `"repartition"` or `"error"`, see
#'   *Details*.
#' @param .diagnose_combination A bool. If `TRUE` (default) the joint group's
#'   single-row (fixed-dose combination) share is computed for
#'   `group_summary$n_single_row_both`.
#'
#' @return A named list:
#' * `ids`: named list of sorted, unique `primaryid` character vectors.
#' * `raw_ids`: report ID sets named by the input definitions, before
#'   splitting or overlap handling.
#' * `membership`: [data.table][data.table::data.table] with one row per scope
#'   report and one logical column per group.
#' * `summary`: [data.table][data.table::data.table] with `group`, `n_reports`,
#'   `share_of_scope`, `n_removed_overlap`, `is_empty`, the definition columns
#'   `drugs`, `roles`, `logic`, `exclude_drugs`, `drug_col`, `n_patterns`, and
#'   `n_single_row_both`. The `source_definitions` list column contains the
#'   named exposure definitions used to build the group, with resolved drug
#'   columns. For an intersection, `drugs` and `n_patterns` cover both
#'   definitions, `roles` are keyed by source group, and `logic = "and"`
#'   describes the intersection of the two definitions. The scalar
#'   `drug_col` is `NA` if the definitions use different columns.
#'   For a group derived from a single definition,
#'   `n_removed_overlap` counts reports removed from that definition. The
#'   intersection and `Neither` have zero removals.
#' * `overlaps`: [data.table][data.table::data.table] with `group1`, `group2`,
#'   `n_intersection` and `strategy`.
#' * `scope_ids`: report IDs in the study scope.
#'
#' @examples
#' data <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
#' data <- faers_dedup(data)
#'
#' # three overlapping exposure groups compared against the same scope
#' groups <- faers_group_ids(
#'     data,
#'     list(
#'         insulin = "insulin|humalog",
#'         biguanide = "metformin",
#'         sulfonylurea = "glyburide|glipizide"
#'     )
#' )
#' groups$summary
#' groups$overlaps
#'
#' # two drugs crossed into A-only / B-only / A+B / Neither
#' joint <- faers_group_ids(
#'     data, list(A = "insulin|humalog", B = "metformin"),
#'     .split = "joint"
#' )
#' lengths(joint$ids)
#'
#' @seealso [faers_compare()] for contingency tables, signal metrics and
#'   matrices built on these groups; [faers_reports_with()] for the exposure
#'   matching semantics; [faers_restrict()] for fixing the scope first.
#' @export
faers_group_ids <- function(.object,
                            .groups,
                            .split = c("none", "joint"),
                            .overlap = c("keep", "exclude", "repartition", "error"),
                            .diagnose_combination = TRUE) {
    .split <- match.arg(.split)
    .overlap <- match.arg(.overlap)
    groups <- .compare_build_groups(.object, .groups, .split, .overlap,
        .diagnose_combination
    )
    groups$raw_ids_by_group <- NULL
    groups
}

.compare_build_groups <- function(.object, .groups, .split, .overlap,
                                  .diagnose_combination) {
    assert_bool(.diagnose_combination)
    assert_s3_class(.object, "FAERSascii")
    if (!.object@standardization) {
        cli::cli_abort("{.arg .object} must be standardized using {.fn faers_standardize}")
    }

    defs <- .compare_check_groups(.groups, .split)
    scope_ids <- sort(unique(as.character(faers_primaryid(.object))))
    n_scope <- length(scope_ids)

    # ---- (1) raw exposure sets, one per definition ----
    raw_defs <- lapply(defs, function(def) {
        as.character(.compare_group_ids(.object, def))
    })
    names(raw_defs) <- names(defs)

    # ---- (2) joint split, or overlap reconciliation ----
    # Record source definitions in the same order as the resulting sets.
    # Group names are labels and can contain the separators used in a split.
    if (identical(.split, "joint")) {
        ids <- .compare_split_joint(raw_defs, scope_ids)
        sources <- list(1L, 2L, c(1L, 2L), integer())
    } else {
        ids <- .compare_apply_overlap(raw_defs, .overlap)
        sources <- if (identical(.overlap, "repartition")) {
            list(1L, 2L, c(1L, 2L))
        } else {
            as.list(seq_along(defs))
        }
    }

    # Compare each final group to its source definition. The intersection and
    # Neither are new sets, so their baseline is their own membership.
    raw_ids_by_group <- stats::setNames(lapply(seq_along(ids), function(i) {
        source <- sources[[i]]
        if (length(source) == 1L) {
            raw_defs[[source]]
        } else {
            ids[[i]]
        }
    }), names(ids))

    # ---- (3) pairwise overlap diagnosis (on the definition-level raw sets) ----
    overlaps <- .compare_overlap_table(raw_defs, .overlap)

    # ---- (4) single-row (fixed-dose) share of the joint group ----
    single_row <- 0L
    if (isTRUE(.diagnose_combination) && length(defs) == 2L) {
        single_row <- .compare_single_row_both(
            .object, defs, intersect(raw_defs[[1L]], raw_defs[[2L]])
        )
    }

    # ---- (5) summary ----
    summary_defs <- lapply(defs, function(def) {
        def$.drug_col <- .compare_group_drug_col(.object, def)
        def[.COMPARE_GROUP_KEYS]
    })
    summary <- data.table::rbindlist(lapply(seq_along(ids), function(i) {
        source <- summary_defs[sources[[i]]]
        n <- length(ids[[i]])
        is_both <- length(source) == 2L
        drugs <- as.character(unlist(lapply(source, `[[`, ".drugs"), use.names = FALSE))
        exclusions <- as.character(unlist(lapply(source, `[[`, ".exclude_drugs"), use.names = FALSE))
        roles <- lapply(source, `[[`, ".roles")
        if (length(roles) == 1L) {
            roles <- roles[[1L]]
        }
        drug_cols <- unique(vapply(source, `[[`, character(1L), ".drug_col"))
        data.table::data.table(
            group = names(ids)[[i]],
            n_reports = as.integer(n),
            share_of_scope = if (n_scope) n / n_scope else NA_real_,
            n_removed_overlap = as.integer(length(raw_ids_by_group[[i]]) - n),
            is_empty = n == 0L,
            drugs = list(drugs),
            roles = list(roles),
            logic = if (is_both) "and" else if (length(source)) source[[1L]]$.logic else NA_character_,
            exclude_drugs = list(unique(exclusions)),
            drug_col = if (length(drug_cols) == 1L) drug_cols[[1L]] else NA_character_,
            n_patterns = length(drugs),
            n_single_row_both = as.integer(if (is_both) single_row else 0L),
            source_definitions = list(source)
        )
    }))

    # ---- (6) membership matrix ----
    membership <- data.table::data.table(primaryid = scope_ids)
    for (nm in names(ids)) {
        membership[[nm]] <- scope_ids %chin% ids[[nm]]
    }

    list(
        ids = ids,
        raw_ids = raw_defs,
        membership = membership,
        summary = summary[],
        overlaps = overlaps[],
        scope_ids = scope_ids,
        raw_ids_by_group = raw_ids_by_group
    )
}
