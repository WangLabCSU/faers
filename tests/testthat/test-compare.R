test_that("background mode conserves the scope in every contingency table", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data,
        list(A = "insulin|humalog", B = "metformin", C = "glyburide"),
        .events = "pt"
    ))

    sr <- res$signal_results
    expect_true(all(sr$a >= 0L & sr$b >= 0L & sr$c >= 0L & sr$d >= 0L))
    expect_true(all(sr$a + sr$b + sr$c + sr$d == sr$n_total))
    # Every background comparison uses the same scope.
    expect_identical(res$metadata$scope$n_reports_universe, 200L)
    expect_true(all(sr$n_total == 200L))
    expect_true(all(sr$mode == "background"))
    expect_true(all(sr$background_id == "rest_of_scope"))
    # The report count of a group bounds its event counts.
    expect_true(all(sr$a <= sr$n_total))

    # One comparison per group.
    expect_identical(
        res$comparisons$comparison_id,
        paste0(c("A", "B", "C"), "_vs_rest_of_scope")
    )
    expect_identical(res$comparisons$n1, as.integer(res$group_summary$n_reports))
    expect_identical(res$comparisons$n0, 200L - as.integer(res$group_summary$n_reports))
})

test_that("the event profile counts reports once, with an explicit denominator", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"), .events = "pt"
    ))
    prof <- res$event_profile

    expect_true(all(prof$N >= 0L))
    expect_true(all(prof$N <= prof$denominator))
    expect_true(all(prof$denominator == res$group_summary$n_reports[
        match(prof$group, res$group_summary$group)
    ]))
    expect_equal(prof$share, prof$N / prof$denominator)
    expect_true(all(prof$universe_N >= prof$N))
    expect_equal(prof$universe_share, prof$universe_N / 200)

    # Every group is laid out over the same candidate events, so the profile is
    # a complete grid that lines up with `signal_results`.
    expect_identical(
        sort(unique(prof$pt)),
        sort(unique(res$signal_results$pt))
    )
    expect_identical(
        prof[, .N, by = "group"]$N,
        rep.int(length(unique(prof$pt)), length(unique(prof$group)))
    )

    # Duplicated drug/event rows must not inflate a count: the group count has
    # one row per (report, event) pair, and `faers_counts()` counts DISTINCT
    # primaryid, so a report with two matching drug rows still counts once.
    # The profile is the full grid (`N = 0` for absent events), so only the rows
    # the group actually has are compared against the sparse count table.
    ids_a <- compare_quiet(faers_reports_with(data, "insulin|humalog", .return = "ids"))
    direct <- compare_quiet(faers_counts(faers_keep(data, ids_a), .events = "pt"))
    sub <- prof[group == "A" & N > 0L][order(pt)]
    expect_identical(sub$N, direct[order(pt)]$N[match(sub$pt, direct[order(pt)]$pt)])
    expect_identical(nrow(sub), nrow(direct))
    expect_true(all(prof[group == "A" & !pt %chin% direct$pt]$N == 0L))
})

test_that("a batch result matches the equivalent single call byte for byte", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data,
        list(insulin = "insulin|humalog", metformin = "metformin"),
        .events = "pt", .methods = c("ror", "prr")
    ))

    ids <- compare_quiet(faers_reports_with(data, "insulin|humalog", .return = "ids"))
    single <- compare_quiet(faers_phv_signal(
        faers_keep(data, ids),
        .events = "pt", .full = data, .methods = c("ror", "prr")
    ))

    batch <- res$signal_results[
        comparison_id == "insulin_vs_rest_of_scope",
        names(single),
        with = FALSE
    ]
    # The only extra column `faers_compare()` adds is `status`, and the
    # `n_below_min_a` attribute tracks a filter that did not drop any row here.
    data.table::setattr(batch, "n_below_min_a", NULL)
    data.table::setkey(batch, pt)
    data.table::setkey(single, pt)
    expect_identical(batch, single)
    expect_identical(res$metadata$event_scope$n_below_min_a, 0L)

    # The cached-scope identity the implementation relies on: comparing a group
    # against the whole scope equals comparing it against the scope minus itself.
    full <- compare_quiet(faers_phv_table(faers_keep(data, ids), .events = "pt", .full = data))
    rest <- compare_quiet(faers_phv_table(
        faers_keep(data, ids),
        .events = "pt", .object2 = faers_keep(data, ids, invert = TRUE)
    ))
    expect_identical(full, rest)
})

test_that("pairwise mode uses the two selected groups as its background", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data,
        list(A = "insulin|humalog", B = "metformin"),
        .mode = "pairwise", .comparisons = list(c("A", "B")),
        .overlap = "exclude", .events = "pt"
    ))

    expect_identical(nrow(res$comparisons), 1L)
    expect_identical(res$comparisons$mode, "pairwise")
    expect_identical(res$comparisons$group1, "A")
    expect_identical(res$comparisons$group2, "B")
    expect_identical(res$comparisons$background_id, "B")
    expect_identical(res$comparisons$comparison_id, "A_vs_B")

    sr <- res$signal_results
    expect_true(all(sr$mode == "pairwise"))
    # Pairwise tables conserve |G1| + |G2|, not the scope.
    expect_true(all(sr$a + sr$b + sr$c + sr$d == sr$n_total))
    expect_identical(
        unique(sr$n_total), res$comparisons$n1 + res$comparisons$n0
    )
    expect_true(all(sr$a <= res$comparisons$n1))
    expect_true(all(sr$c <= res$comparisons$n0))

    # Exclusion is applied before the table is built, and the removal is recorded.
    groups <- compare_quiet(faers_group_ids(
        data, list(A = "insulin|humalog", B = "metformin"), .overlap = "exclude"
    ))
    expect_identical(res$comparisons$n1, as.integer(length(groups$ids$A)))
    expect_identical(res$comparisons$n0, as.integer(length(groups$ids$B)))
    # Two reports matched both drugs; "exclude" took them out of both arms.
    raw <- compare_quiet(faers_group_ids(data, list(A = "insulin|humalog", B = "metformin")))
    shared <- length(intersect(raw$ids$A, raw$ids$B))
    expect_identical(shared, 2L)
    expect_identical(res$comparisons$n_overlap_removed, as.integer(shared))
    expect_identical(res$comparisons$n_overlap_shared, as.integer(shared))
    # The tables are built on the reconciled, disjoint groups.
    expect_identical(length(intersect(groups$ids$A, groups$ids$B)), 0L)
})

test_that("pairwise mode refuses overlapping groups unless reconciled", {
    data <- compare_std_data()
    expect_error(
        compare_quiet(faers_compare(
            data,
            list(A = "insulin|humalog", B = "metformin"),
            .mode = "pairwise", .comparisons = list(c("A", "B")),
            .events = "pt"
        )),
        "disjoint|overlap"
    )
    expect_error(
        compare_quiet(faers_compare(
            data,
            list(A = "insulin|humalog", B = "metformin"),
            .mode = "pairwise", .comparisons = list(c("A", "B")),
            .overlap = "error", .events = "pt"
        ))
    )
})

test_that("pairwise mode never runs implicitly", {
    data <- compare_std_data()
    # Every pair changes the denominator, so an implied batch would invite
    # uncorrected multiple comparisons; the plan must be stated explicitly.
    expect_error(
        compare_quiet(faers_compare(
            data, list(A = "insulin|humalog", B = "metformin"),
            .mode = "pairwise", .events = "pt"
        )),
        "comparison"
    )
})

test_that("joint groups are mutually exclusive and give a four-way partition", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data,
        list(A = "insulin|humalog", B = "metformin"),
        .split = "joint", .events = "pt"
    ))

    expect_identical(
        res$group_summary$group, c("A-only", "B-only", "A+B", "Neither")
    )
    expect_identical(sum(res$group_summary$n_reports), 200L)
    # Background mode still uses the scope as the common denominator.
    sr <- res$signal_results
    expect_true(all(sr$a + sr$b + sr$c + sr$d == 200L))
    # A joint group can be contrasted with a single-only group.
    expect_identical(
        sort(unique(res$group_summary$group)),
        sort(c("A-only", "B-only", "A+B", "Neither"))
    )
})

test_that("status flags zero cells and empty groups without altering the numbers", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"), .events = "pt"
    ))
    sr <- res$signal_results

    expect_true(all(sr$status %chin% c("ok", "zero_cell", "empty_group")))
    expect_true(any(sr$status == "ok"))
    # A zero cell is flagged but the metric is left untouched: `a = 0` yields
    # `ror = 0` and `ror_ci_high = NaN`, exactly as `phv_ror()` produces.
    zero <- sr[a == 0L]
    if (nrow(zero)) {
        expect_true(all(zero$ror == 0))
        expect_true(all(is.nan(zero$ror_ci_high)))
        expect_true(all(zero$status == "zero_cell"))
    }
    # No cell is imputed: the raw zero is still visible in the table.
    expect_true(any(sr$a == 0L))
})

test_that("an empty group is legal and flagged rather than an error", {
    data <- compare_std_data()
    expect_no_error(
        res <- compare_quiet(faers_compare(
            data,
            list(A = "insulin|humalog", empty = "no_such_drug_zzz"),
            .events = "pt"
        ))
    )
    res <- compare_quiet(faers_compare(
        data,
        list(A = "insulin|humalog", empty = "no_such_drug_zzz"),
        .events = "pt"
    ))
    expect_true(res$group_summary[group == "empty", is_empty])
    sr <- res$signal_results[comparison_id == "empty_vs_rest_of_scope"]
    expect_true(all(sr$a == 0L))
    expect_true(all(sr$b == 0L))
    expect_true(all(sr$status == "empty_group"))
})

test_that(".min_a filters rows and reports how many were dropped", {
    data <- compare_std_data()
    loose <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"), .events = "pt"
    ))
    strict <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"),
        .events = "pt", .min_a = 2L
    ))

    expect_true(nrow(strict$signal_results) <= nrow(loose$signal_results))
    expect_true(all(strict$signal_results$a >= 2L))
    expect_identical(strict$metadata$event_scope$min_a, 2L)
    expect_identical(
        strict$metadata$event_scope$n_below_min_a,
        as.integer(sum(loose$signal_results$a < 2L))
    )
})

test_that("matrices are event by group, with a signal matrix per comparison", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data,
        list(A = "insulin|humalog", B = "metformin", C = "glyburide"),
        .events = "pt"
    ))
    m <- res$matrices

    expect_named(m, c("count", "share", "signal"))
    expect_identical(names(m$count), c("pt", "A", "B", "C"))
    expect_identical(names(m$share), c("pt", "A", "B", "C"))
    # The signal matrix is keyed on comparisons, not groups.
    expect_identical(
        names(m$signal),
        c("pt", res$comparisons$comparison_id)
    )
    # Counts are integers and shares are proportions of the group denominator.
    expect_true(all(vapply(m$count[, -1L], is.integer, logical(1L))))
    expect_true(all(unlist(m$share[, -1L]) >= 0 & unlist(m$share[, -1L]) <= 1))
    # The count matrix agrees with the profile (which is a complete grid).
    for (g in c("A", "B", "C")) {
        expect_identical(
            m$count[[g]],
            res$event_profile[group == g]$N[match(m$count$pt, res$event_profile[group == g]$pt)]
        )
    }
})

test_that(".matrix_metric selects the signal matrix column", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"),
        .events = "pt", .matrix_metric = "ror"
    ))
    expect_true("ror" %chin% names(res$matrices))
    expect_identical(res$metadata$matrix_metric, "ror")

    expect_error(
        compare_quiet(faers_compare(
            data, list(A = "insulin|humalog", B = "metformin"),
            .events = "pt", .matrix_metric = "not_a_metric"
        )),
        "not_a_metric"
    )
})

test_that(".top_n truncates the display layer only and flags it", {
    data <- compare_std_data()
    full <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"), .events = "pt"
    ))
    trunc <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"),
        .events = "pt", .top_n = 5L
    ))

    # `signal_results` is never truncated.
    expect_identical(nrow(trunc$signal_results), nrow(full$signal_results))
    # The profile keeps the union of each group's top-5 events and stays a
    # complete grid, so every group shows the same (unioned) event set.
    expect_true(nrow(trunc$event_profile) < nrow(full$event_profile))
    n_groups <- length(unique(trunc$event_profile$group))
    expect_true(length(unique(trunc$event_profile$pt)) <= 5L * n_groups)
    expect_identical(
        trunc$event_profile[, .N, by = "group"]$N,
        rep.int(length(unique(trunc$event_profile$pt)), n_groups)
    )
    # Every retained event is a top-5 event of at least one group.
    top5 <- full$event_profile[
        order(-N, pt), utils::head(.SD, 5L), by = "group"
    ]$pt
    expect_true(all(unique(trunc$event_profile$pt) %chin% top5))
    # Flagged in both the attribute and the metadata.
    expect_true(isTRUE(attr(trunc$event_profile, "display_subset")))
    expect_true(trunc$metadata$display_subset$event_profile)
    expect_true(trunc$metadata$display_subset$matrices)
    expect_false(trunc$metadata$display_subset$signal_results)
    expect_identical(trunc$metadata$truncation$top_n, 5L)
    # Untruncated runs are not flagged.
    expect_false(isTRUE(attr(full$event_profile, "display_subset")))
})

test_that(".event_set narrows candidate events without collapsing them", {
    data <- compare_std_data()
    selected <- c("ASTHENIA", "ANAEMIA", "DRUG INTERACTION")
    res <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"),
        .events = "pt", .event_set = selected
    ))

    expect_setequal(unique(res$signal_results$pt), selected)
    expect_setequal(unique(res$event_profile$pt), selected)
    expect_true(all(res$event_profile$pt %chin% selected))
    # One row per event per group, NOT a single composite row.
    expect_identical(nrow(res$signal_results), 2L * length(selected))
    expect_identical(
        nrow(res$event_profile), 2L * length(selected)
    )
    expect_false(any(grepl("Composite", res$signal_results$pt)))
    expect_identical(res$metadata$event_scope$n_selected_events, 3L)

    # A function selector behaves like the character selector.
    by_function <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"),
        .events = "pt", .event_set = function(x) x %in% selected
    ))
    expect_setequal(unique(by_function$signal_results$pt), selected)

    # An empty selection yields well-formed empty tables instead of an error.
    none <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"),
        .events = "pt", .event_set = "NO_SUCH_EVENT_XYZ"
    ))
    expect_identical(nrow(none$signal_results), 0L)
    expect_identical(nrow(none$event_profile), 0L)
    expect_true(all(
        c("group", "pt", "N", "denominator", "share", "universe_N", "universe_share") %chin%
            names(none$event_profile)
    ))
    # Regression: an empty result must carry the same columns and column types
    # as a populated one, metric columns included. It used to be built from a
    # hardcoded name list that mismatched `phv_signal()` (e.g. `chisq` without
    # `chisq_pvalue`), so an explicit `.matrix_metric = "ror"` aborted on an
    # empty selection even though it worked on a populated one.
    expect_identical(names(none$signal_results), names(res$signal_results))
    expect_identical(
        vapply(none$signal_results, function(x) class(x)[[1L]], character(1L)),
        vapply(res$signal_results, function(x) class(x)[[1L]], character(1L))
    )
    # `.min_a` accounting survives the empty branch.
    expect_identical(attr(none$signal_results, "n_below_min_a"), 0L)
    for (metric in c("ror", "expected", "prr")) {
        expect_no_error(compare_quiet(faers_compare(
            data, list(A = "insulin|humalog", B = "metformin"),
            .events = "pt", .event_set = "NO_SUCH_EVENT_XYZ",
            .matrix_metric = metric
        )))
    }
    # A metric that really is unknown is still rejected, empty selection or not.
    expect_error(
        compare_quiet(faers_compare(
            data, list(A = "insulin|humalog", B = "metformin"),
            .events = "pt", .event_set = "NO_SUCH_EVENT_XYZ",
            .matrix_metric = "no_such_metric"
        )),
        "is not a column"
    )

    # A selector of the wrong shape is rejected.
    expect_error(
        compare_quiet(faers_compare(
            data, list(A = "insulin|humalog", B = "metformin"),
            .events = "pt", .event_set = c(TRUE, FALSE)
        )),
        "one value per candidate event"
    )
    expect_error(
        compare_quiet(faers_compare(
            data, list(A = "insulin|humalog", B = "metformin"),
            .events = "pt", .event_set = 1:3
        )),
        "character vector, a logical vector, or a function"
    )
})

test_that("multi-column events use a collapsed event_label in the matrices", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"),
        .events = c("pt", "soc_name"), .min_a = 1L
    ))

    expect_true("event_label" %chin% names(res$event_profile))
    expect_identical(names(res$matrices$count)[[1L]], "event_label")
    expect_identical(names(res$matrices$signal)[[1L]], "event_label")
    # The label is the documented paste of the event columns.
    expect_true(all(grepl(" | ", res$event_profile$event_label, fixed = TRUE)))
    # It also survives truncation.
    trunc <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"),
        .events = c("pt", "soc_name"), .min_a = 1L, .top_n = 4L
    ))
    expect_identical(names(trunc$matrices$count)[[1L]], "event_label")
})

test_that("metadata documents the scope, groups and algorithms", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"), .events = "pt"
    ))
    md <- res$metadata

    expect_true(all(c(
        "call", "faers_version", "r_version", "timestamp", "mode", "split",
        "overlap_strategy", "groups", "scope", "events", "event_scope",
        "methods", "denominator_definition", "background_definition",
        "comparisons", "missing", "truncation", "display_subset",
        "interaction_test"
    ) %chin% names(md)))
    expect_identical(md$mode, "background")
    expect_identical(md$split, "none")
    expect_identical(md$overlap_strategy, "keep")
    expect_identical(md$scope$n_reports_universe, 200L)
    expect_true(md$scope$standardization)
    expect_true(md$scope$deduplication)
    expect_identical(md$events, "pt")
    expect_identical(md$methods, c("ror", "prr"))
    expect_true(nzchar(md$faers_version))
    expect_true(nzchar(md$r_version))
    expect_true(nzchar(md$denominator_definition))
    expect_true(nzchar(md$background_definition))
    # The joint-group result is explicitly not a formal interaction test.
    expect_match(md$interaction_test, "interaction")
    # `.groups` is echoed verbatim for reproducibility.
    expect_identical(md$groups, list(A = "insulin|humalog", B = "metformin"))
})

test_that("faers_compare() requires standardized, de-duplicated input", {
    raw <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
    expect_error(
        compare_quiet(faers_compare(raw, list(A = "insulin"), .events = "pt")),
        "de-duplicated"
    )
})

test_that("faers_compare() validates its arguments", {
    data <- compare_std_data()
    groups <- list(A = "insulin|humalog")

    expect_error(compare_quiet(faers_compare(data, groups, .events = character())))
    expect_error(compare_quiet(faers_compare(data, groups, .events = "not_a_column")))
    expect_error(compare_quiet(faers_compare(data, groups, .min_a = -1L)))
    expect_error(compare_quiet(faers_compare(data, groups, .methods = "nope")))
    expect_error(compare_quiet(faers_compare(data, groups, .top_n = 0L)))
    expect_error(compare_quiet(faers_compare(data, groups, .matrix_metric = "")))
})

test_that("faers_compare() is identical on both backends", {
    skip_if_not_installed("duckdb")
    skip_if_not_installed("DBI")

    data <- compare_std_data()
    db <- compare_db_twin(data)
    on.exit(DBI::dbDisconnect(db@db@con, shutdown = TRUE), add = TRUE)

    groups <- list(A = "insulin|humalog", B = "metformin", C = "glyburide")
    cases <- list(
        list(defs = groups, args = list(.events = "pt")),
        list(defs = groups, args = list(.events = "soc_name")),
        list(defs = groups, args = list(.events = "pt", .event_set = c("ASTHENIA", "ANAEMIA"))),
        list(defs = groups, args = list(.events = "pt", .top_n = 4L, .min_a = 1L)),
        list(defs = groups, args = list(.events = c("pt", "soc_name"), .min_a = 1L)),
        list(
            defs = groups[1:2],
            args = list(
                .split = "joint", .mode = "pairwise",
                .comparisons = list(c("A-only", "B-only"), c("A+B", "A-only")),
                .events = "pt"
            )
        )
    )
    for (case in cases) {
        mem <- compare_quiet(do.call(
            faers_compare, c(list(data, case$defs), case$args)
        ))
        dbo <- compare_quiet(do.call(
            faers_compare, c(list(db, case$defs), case$args)
        ))
        for (part in c("group_summary", "event_profile", "comparisons", "group_overlaps")) {
            expect_identical(mem[[part]], dbo[[part]])
        }
        expect_identical(
            mem$signal_results[, .(comparison_id, a, b, c, d, n_total)],
            dbo$signal_results[, .(comparison_id, a, b, c, d, n_total)]
        )
        expect_identical(mem$matrices$count, dbo$matrices$count)
        expect_identical(mem$matrices$share, dbo$matrices$share)
    }
})
