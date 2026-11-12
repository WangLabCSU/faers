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

test_that("derived comparisons count shared reports absent from both arms", {
    data <- compare_std_data()
    groups <- list("A & partner" = "insulin|humalog", "B+other" = "metformin")
    raw <- faers_group_ids(data, groups)$raw_ids
    n_shared <- length(intersect(raw[[1L]], raw[[2L]]))
    expect_gt(n_shared, 0L)

    for (split in c("none", "joint")) {
        overlap <- if (split == "joint") "keep" else "repartition"
        both <- if (split == "joint") "A & partner+B+other" else "A & partner&B+other"
        only <- paste0(names(groups), "-only")
        grouped <- faers_group_ids(data, groups, .split = split, .overlap = overlap)
        expect_identical(grouped$raw_ids, raw)
        expect_identical(
            grouped$summary$n_removed_overlap[match(c(only, both), grouped$summary$group)],
            c(n_shared, n_shared, 0L)
        )
        expect_identical(grouped$summary$drugs[[1L]], groups[[1L]])
        expect_identical(grouped$summary$drugs[[2L]], groups[[2L]])

        res <- compare_quiet(faers_compare(
            data, groups, .split = split, .overlap = overlap,
            .mode = "pairwise", .comparisons = list(only, c(both, only[[1L]]))
        ))
        expect_identical(res$comparisons$n_overlap_shared, rep(n_shared, 2L))
        expect_identical(res$comparisons$n_overlap_removed, c(n_shared, 0L))
    }
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

test_that("an empty side is flagged in either comparison direction", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", empty = "no_such_drug_zzz"),
        .mode = "pairwise",
        .comparisons = list(c("A", "empty"), c("empty", "A"))
    ))
    sr <- res$signal_results
    expect_setequal(unique(sr$comparison_id), c("A_vs_empty", "empty_vs_A"))
    expect_true(all(sr$status == "empty_group"))
    expect_true(all(sr[comparison_id == "A_vs_empty", c + d] == 0L))
    expect_true(all(sr[comparison_id == "empty_vs_A", a + b] == 0L))
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

test_that("repeated or colliding comparisons are refused, never aggregated", {
    data <- compare_std_data()
    groups <- list(A = "insulin|humalog", B = "metformin", C = "glyburide")

    # The same pair twice used to yield one shared `comparison_id`, so the cast
    # that builds `matrices$signal` aggregated the metric into a row count and
    # returned a plausible-looking matrix of 2s.
    expect_error(
        compare_quiet(faers_compare(data, groups,
            .mode = "pairwise", .overlap = "exclude",
            .comparisons = list(c("A", "B"), c("A", "B")), .events = "pt"
        )),
        "must define distinct comparisons"
    )
    # Two *different* pairs can also collapse onto one id when a group name
    # already contains the separator: "A_vs_B" vs "A" and "A" vs "B_vs_A".
    named <- list(A = "insulin|humalog", B = "metformin",
        "A_vs_B" = "glyburide", "B_vs_A" = "glipizide")
    expect_error(
        compare_quiet(faers_compare(data, named,
            .mode = "pairwise", .overlap = "exclude",
            .comparisons = list(c("A_vs_B", "A"), c("A", "B_vs_A")), .events = "pt"
        )),
        "A_vs_B_vs_A"
    )
    # Order matters, so the reversed pair is a genuinely different comparison.
    ok <- compare_quiet(faers_compare(data, groups,
        .mode = "pairwise", .overlap = "exclude",
        .comparisons = list(c("A", "B"), c("B", "A")), .events = "pt"
    ))
    expect_identical(
        setdiff(names(ok$matrices$signal), "pt"),
        c("A_vs_B", "B_vs_A")
    )
    expect_length(unique(ok$comparisons$comparison_id), 2L)
})

test_that("the matrix cast refuses a repeated key instead of aggregating", {
    # Defence in depth: even if a caller ever hands `dcast()` a repeated key,
    # the matrix must fail loudly rather than turn values into row counts.
    bad <- data.table::data.table(
        event_label = c("x", "x"), comparison_id = c("c", "c"), value = c(1, 2)
    )
    expect_error(
        faers:::.compare_matrix(bad, "event_label", "comparison_id", "value", "c"),
        "is not unique in the matrix source"
    )
    # A zero-row source is legal (an empty group) and still casts.
    empty <- bad[0L]
    expect_identical(
        nrow(faers:::.compare_matrix(empty, "event_label", "comparison_id", "value", "c")),
        0L
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

test_that("display matrices retain comparisons without any signal rows", {
    data <- compare_std_data()
    groups <- list(A = "insulin|humalog", B = "no_matching_drug_xyz")

    some <- compare_quiet(faers_compare(
        data, groups, .top_n = 2L, .min_a = 1L, .matrix_metric = "expected"
    ))
    ids <- some$comparisons$comparison_id
    signal <- some$matrices$expected
    expect_named(signal, c("pt", ids))
    expect_identical(signal$pt, some$matrices$count$pt)
    expect_identical(signal$pt, some$matrices$share$pt)
    expect_false(ids[[2L]] %in% some$signal_results$comparison_id)
    expect_identical(signal[[ids[[2L]]]], rep(NA_real_, nrow(signal)))
    expect_true(any(is.finite(signal[[ids[[1L]]]])))

    none <- compare_quiet(faers_compare(
        data, groups, .top_n = 2L, .min_a = length(faers_primaryid(data)) + 1L
    ))
    expect_equal(nrow(none$signal_results), 0L)
    expect_gt(nrow(none$matrices$count), 0L)
    expect_named(none$matrices$signal, c("pt", ids))
    expect_identical(none$matrices$signal$pt, none$matrices$count$pt)
    for (id in ids) {
        expect_identical(none$matrices$signal[[id]], rep(NA_real_, nrow(none$matrices$count)))
    }
})

test_that(".top_n aligns matrices for numeric event columns", {
    data <- compare_std_data()
    res <- compare_quiet(faers_compare(
        data, list(A = "insulin|humalog", B = "metformin"),
        .events = "age_in_years", .field = "demo", .top_n = 2L
    ))
    labels <- res$matrices$count$age_in_years
    expect_true(is.numeric(labels))
    expect_gt(length(labels), 0L)
    expect_identical(res$matrices$share$age_in_years, labels)
    expect_identical(res$matrices$signal$age_in_years, labels)
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

test_that(".event_set selects numeric and missing event values", {
    data <- compare_std_data()
    groups <- list(A = "insulin|humalog", B = "metformin")
    demo <- data.table::copy(faers_get(data, "demo"))
    demo[1L, age_in_years := NA_real_]
    data@data$demo <- demo
    ages <- demo$age_in_years
    compare_age <- function(object, selector, na_rm = FALSE) {
        compare_quiet(faers_compare(
            object, groups, .field = "demo", .events = "age_in_years",
            .event_set = selector, .na.rm = na_rm
        ))
    }

    older <- compare_age(data, function(x) x >= 50)
    expect_setequal(older$event_profile$age_in_years, ages[!is.na(ages) & ages >= 50])
    expect_identical(older$metadata$event_scope$n_candidate_events, length(unique(ages)))

    age <- ages[!is.na(ages)][[1L]]
    expect_identical(unique(compare_age(data, as.character(age))$signal_results$age_in_years), age)
    missing_age <- compare_age(data, is.na)
    expect_identical(unique(missing_age$event_profile$age_in_years), NA_real_)
    expect_identical(missing_age$metadata$event_scope$n_selected_events, 1L)
    expect_identical(
        unique(missing_age$event_profile$universe_N),
        sum(is.na(ages))
    )
    dropped <- compare_age(data, is.na, TRUE)
    expect_identical(dropped$event_profile$age_in_years, numeric())
    expect_identical(dropped$signal_results$age_in_years, numeric())

    skip_if_not_installed("duckdb")
    skip_if_not_installed("DBI")
    db <- compare_db_twin(data)
    on.exit(DBI::dbDisconnect(db@db@con, shutdown = TRUE), add = TRUE)
    expect_identical(compare_age(db, function(x) x >= 50)$signal_results, older$signal_results)
    expect_identical(compare_age(db, is.na)$event_profile, missing_age$event_profile)
})

test_that(".fn runs once on the study scope and all groups use its rows", {
    data <- compare_std_data()
    groups <- list(A = "insulin|humalog", B = "metformin")
    demo <- faers_get(data, "demo")
    original_sex <- data.table::copy(demo$sex)
    first_id <- sort(demo$primaryid)[[1L]]
    ids <- faers_group_ids(data, groups)$ids
    calls <- 0L
    keep_first <- function(rows, limit) {
        calls <<- calls + 1L
        data.table::set(rows, j = "sex", value = rep("chosen", nrow(rows)))
        rows[order(primaryid)][seq_len(limit)]
    }
    compare_first <- function(object) {
        compare_quiet(faers_compare(
            object, groups, .field = "demo", .events = "sex",
            .fn = keep_first, limit = 1L
        ))
    }
    mem <- compare_first(data)
    expect_identical(calls, 1L)
    expect_identical(faers_get(data, "demo")$sex, original_sex)
    expect_identical(mem$signal_results$sex, rep("chosen", length(groups)))
    expect_identical(mem$signal_results$a, as.integer(vapply(ids, function(x) first_id %in% x, logical(1L))))
    expect_true(all(with(mem$signal_results, a + c == 1L)))
    expect_true(all(with(mem$signal_results, a + b + c + d == n_total)))

    expect_error(compare_quiet(faers_compare(
        data, groups, .fn = function(rows) rows[, .(pt)]
    )), "retain the required columns.*primaryid")
    expect_error(compare_quiet(faers_compare(
        data, groups, .fn = function(rows) rows[1L][, primaryid := "outside-scope"]
    )), "report IDs from the study scope")

    skip_if_not_installed("duckdb")
    skip_if_not_installed("DBI")
    db <- compare_db_twin(data)
    on.exit(DBI::dbDisconnect(db@db@con, shutdown = TRUE), add = TRUE)
    db_res <- compare_first(db)
    expect_identical(calls, 2L)
    expect_identical(db_res$event_profile, mem$event_profile)
    expect_identical(db_res$signal_results, mem$signal_results)
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

test_that("a non-reac .field counts its own event columns", {
    data <- compare_std_data()
    groups <- list(A = "insulin|humalog", B = "metformin")

    # Regression: `.events` was validated against the `reac` allow-list no
    # matter what `.field` said, so every other field aborted before doing any
    # work (with an unrelated cli pluralization error). The validation now
    # follows `.field`, which makes the documented `...` pass-through usable.
    res <- compare_quiet(faers_compare(
        data, groups, .events = "drugname", .field = "drug", .min_a = 0L
    ))
    expect_identical(res$metadata$event_scope$field, "drug")
    expect_identical(res$metadata$event_scope$allowed_columns,
        names(faers_get(data, "drug")))
    expect_true("drugname" %chin% names(res$event_profile))
    expect_setequal(unique(res$event_profile$group), c("A", "B"))

    # The counts still come from `faers_counts()`, so the group column has to
    # reproduce a direct call on the same reports.
    ids_a <- compare_quiet(faers_reports_with(data, "insulin|humalog", .return = "ids"))
    direct <- compare_quiet(faers_counts(
        faers_keep(data, ids_a), .events = "drugname", .field = "drug"
    ))
    prof <- res$event_profile[group == "A" & N > 0L][order(drugname)]
    expect_identical(prof$N, direct[order(drugname)]$N[match(prof$drugname, direct[order(drugname)]$drugname)])
    expect_identical(nrow(prof), nrow(direct))

    # Every field the object carries, not just `reac`, is reachable.
    for (case in list(
        list(field = "demo", events = "sex"),
        list(field = "outc", events = "outc_cod"),
        list(field = "rpsr", events = "rpsr_cod"),
        list(field = "indi", events = "indi_pt")
    )) {
        out <- compare_quiet(faers_compare(
            data, groups, .events = case$events, .field = case$field, .min_a = 0L
        ))
        expect_identical(out$metadata$event_scope$field, case$field)
        expect_true(nrow(out$event_profile) > 0L)
    }
})

test_that("event columns are validated against the chosen .field, in cli style", {
    data <- compare_std_data()
    groups <- list(A = "insulin|humalog")

    # A column of another field is rejected, and the message names both the
    # column and the field it was looked for in.
    expect_error(
        compare_quiet(faers_compare(data, groups, .events = "pt", .field = "drug")),
        "Unknown column of the standardized drug data.*pt"
    )
    # Singular/plural are both grammatical (cli quantity).
    expect_error(
        compare_quiet(faers_compare(data, groups, .events = "nope", .field = "reac")),
        "Unknown column of the standardized reac data"
    )
    expect_error(
        compare_quiet(faers_compare(data, groups, .events = c("nope1", "nope2"))),
        "Unknown columns of the standardized reac data"
    )
    # `reac` keeps its own allow-list: the MedDRA hierarchy columns are joined
    # in by standardization and are not physical `reac` columns.
    expect_error(
        compare_quiet(faers_compare(data, groups, .events = "caseid")),
        "Unknown column"
    )
    hier <- compare_quiet(faers_compare(data, groups, .events = "soc_name", .min_a = 0L))
    expect_true("soc_name" %chin% names(hier$event_profile))

    # `.field` itself is checked before anything else runs.
    expect_error(
        compare_quiet(faers_compare(data, groups, .events = "pt", .field = "nope")),
        "must contain only values from"
    )
    expect_error(
        compare_quiet(faers_compare(data, groups, .events = "pt", .field = c("reac", "drug"))),
        "a single string"
    )
})

test_that("missing-event metadata agrees with counts and denominators", {
    data <- compare_std_data()
    groups <- list(A = "insulin|humalog")
    ids <- suppressMessages(faers_reports_with(data, groups$A, .return = "ids"))
    demo <- data.table::copy(faers_get(data, "demo"))
    demo[primaryid %in% ids[[1L]], sex := NA_character_]
    data@data$demo <- demo

    compare_sex <- function(object, na_rm) {
        compare_quiet(faers_compare(
            object, groups, .field = "demo", .events = "sex", .na.rm = na_rm
        ))
    }
    kept <- compare_sex(data, FALSE)
    dropped <- compare_sex(data, TRUE)
    expect_false(kept$metadata$missing$na_rm)
    expect_true(dropped$metadata$missing$na_rm)
    expect_identical(kept$event_profile[is.na(sex)]$N,
        as.integer(sum(demo$primaryid %in% ids & is.na(demo$sex))))
    expect_gt(kept$event_profile[is.na(sex)]$N, 0L)
    expect_false(anyNA(dropped$event_profile$sex))
    expect_false(anyNA(dropped$signal_results$sex))
    expect_identical(dropped$event_profile$denominator,
        rep(length(ids), nrow(dropped$event_profile)))
    expect_identical(dropped$comparisons$n_total, kept$comparisons$n_total)
    expect_true(all(with(dropped$signal_results, a + b + c + d == n_total)))

    skip_if_not_installed("duckdb")
    skip_if_not_installed("DBI")
    db <- compare_db_twin(data)
    on.exit(DBI::dbDisconnect(db@db@con, shutdown = TRUE), add = TRUE)
    for (na_rm in c(FALSE, TRUE)) {
        db_res <- compare_sex(db, na_rm)
        mem_res <- if (na_rm) dropped else kept
        expect_identical(db_res$metadata$missing, mem_res$metadata$missing)
        expect_identical(db_res$event_profile, mem_res$event_profile)
        expect_identical(db_res$signal_results, mem_res$signal_results)
    }
})

test_that("overlap metadata counts only intersecting pairs and unique reports", {
    data <- compare_std_data()
    groups <- list(
        A = "insulin|humalog", B = "insulin|humalog", C = "insulin|humalog",
        empty = "no_such_drug_zzz"
    )
    res <- compare_quiet(faers_compare(data, groups, .overlap = "exclude"))
    n_shared <- faers_reports_with(data, groups$A, .return = "count")
    expect_gt(n_shared, 0L)
    expect_identical(nrow(res$group_overlaps), 6L)
    expect_identical(res$metadata$overlap_handling$n_shared_pairs, 3L)
    expect_identical(res$metadata$overlap_handling$n_reports_shared, n_shared)

    disjoint <- compare_quiet(faers_compare(
        data, groups[c("A", "empty")], .mode = "pairwise",
        .comparisons = list(c("A", "empty"))
    ))
    expect_identical(disjoint$metadata$overlap_handling$n_shared_pairs, 0L)
    expect_identical(disjoint$metadata$overlap_handling$n_reports_shared, 0L)
    expect_true(disjoint$metadata$overlap_handling$pairwise)
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
    expect_error(compare_quiet(faers_compare(data, groups, .field = "pt")))
})

test_that("matrix event columns cannot collide with group or comparison names", {
    data <- compare_std_data()
    expect_error(compare_quiet(faers_compare(data, list(pt = "insulin"), .events = "pt")),
        "group name.*conflicts with the matrix event column")

    demo <- data.table::copy(faers_get(data, "demo"))
    demo[, A_vs_rest_of_scope := sex]
    data@data$demo <- demo
    expect_error(
        compare_quiet(faers_compare(
            data, list(A = "insulin"), .field = "demo", .events = "A_vs_rest_of_scope"
        )),
        "comparison_id.*conflicts with the matrix event column"
    )
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
        list(defs = groups, args = list(.events = "drugname", .field = "drug", .min_a = 0L)),
        list(defs = groups, args = list(.events = "sex", .field = "demo", .min_a = 0L)),
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

test_that("faers_compare() never writes over a user table on duckdb", {
    skip_if_not_installed("duckdb")
    skip_if_not_installed("DBI")

    data <- compare_std_data()
    db <- compare_db_twin(data)
    con <- db@db@con
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

    # The group ids used to be registered under one fixed, guessable name, and
    # `dbWriteTable(temporary = TRUE)` REPLACES a same-named permanent table
    # instead of shadowing it. Both the old fixed name and the exact name the
    # counter hands out first are therefore occupied by real user tables here.
    fixed <- "__faers_compare_ids__"
    first <- sprintf("__faers_compare_ids_%d_1__", Sys.getpid())
    for (nm in c(fixed, first)) {
        DBI::dbExecute(con, sprintf('CREATE TABLE "%s" (sentinel VARCHAR)', nm))
        DBI::dbExecute(con, sprintf('INSERT INTO "%s" VALUES (\'keep me\')', nm))
    }
    before <- DBI::dbGetQuery(
        con, "SELECT table_name FROM duckdb_tables()"
    )$table_name
    n_before <- length(DBI::dbListTables(con))

    groups <- list(A = "insulin|humalog", B = "metformin")
    mem <- compare_quiet(faers_compare(data, groups, .events = "pt"))
    dbo <- compare_quiet(faers_compare(db, groups, .events = "pt"))

    # The comparison still produces the in-memory answer ...
    expect_identical(mem$matrices$count, dbo$matrices$count)
    expect_identical(mem$event_profile, dbo$event_profile)

    # ... and both user tables are untouched, contents included.
    expect_true(all(c(fixed, first) %in% DBI::dbGetQuery(
        con, "SELECT table_name FROM duckdb_tables()"
    )$table_name))
    for (nm in c(fixed, first)) {
        expect_identical(
            DBI::dbGetQuery(con, sprintf('SELECT sentinel FROM "%s"', nm))$sentinel,
            "keep me"
        )
    }
    expect_setequal(
        setdiff(DBI::dbListTables(con), before),
        character()
    )
    expect_length(DBI::dbListTables(con), n_before)
})
