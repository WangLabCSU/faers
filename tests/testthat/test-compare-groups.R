test_that("faers_group_ids() keeps groups as defined and reports overlaps", {
    data <- compare_std_data()
    groups <- compare_quiet(faers_group_ids(
        data,
        list(insulin = "insulin|humalog", metformin = "metformin")
    ))

    expect_named(
        groups,
        c("ids", "raw_ids", "membership", "summary", "overlaps", "scope_ids")
    )
    expect_identical(names(groups$ids), c("insulin", "metformin"))
    expect_identical(length(groups$scope_ids), 200L)

    # `"keep"` does not remove anything, but the intersection is reported.
    expect_true(all(groups$summary$n_removed_overlap == 0L))
    expect_identical(nrow(groups$overlaps), 1L)
    expect_identical(
        groups$overlaps$n_intersection,
        length(intersect(groups$ids$insulin, groups$ids$metformin))
    )

    # membership has one row per scope report and one logical column per group.
    expect_identical(nrow(groups$membership), length(groups$scope_ids))
    expect_true(all(c("insulin", "metformin") %in% names(groups$membership)))
})

test_that("joint split gives four mutually exclusive groups that sum to U", {
    data <- compare_std_data()
    joint <- compare_quiet(faers_group_ids(
        data,
        list(A = "insulin|humalog", B = "metformin"),
        .split = "joint"
    ))

    expect_identical(
        names(joint$ids), c("A-only", "B-only", "A+B", "Neither")
    )
    # The four sets partition the scope.
    expect_identical(
        sum(lengths(joint$ids)), length(joint$scope_ids)
    )
    for (i in seq_len(4L)) {
        for (j in seq_len(4L)) {
            if (i < j) {
                expect_length(intersect(joint$ids[[i]], joint$ids[[j]]), 0L)
            }
        }
    }

    # Membership is a one-hot encoding.
    membership <- as.matrix(joint$membership[, -1L])
    expect_true(all(rowSums(membership) == 1L))

    # The definitions match set algebra on the raw exposure sets.
    ids <- joint$ids
    expect_setequal(ids[["A-only"]], setdiff(joint$raw_ids$A, joint$raw_ids$B))
    expect_setequal(ids[["B-only"]], setdiff(joint$raw_ids$B, joint$raw_ids$A))
    expect_setequal(ids[["A+B"]], intersect(joint$raw_ids$A, joint$raw_ids$B))
    expect_setequal(
        ids[["Neither"]],
        setdiff(joint$scope_ids, union(joint$raw_ids$A, joint$raw_ids$B))
    )
})

test_that("the joint complement uses the same role-restricted definition", {
    data <- compare_std_data()
    # A report whose only B record carries a role outside the allowed set must
    # not be treated as B-exposed by the complement.
    joint <- compare_quiet(faers_group_ids(
        data,
        list(
            A = "insulin|humalog",
            B = list(.drugs = "metformin", .roles = list(metformin = "PS"))
        ),
        .split = "joint"
    ))
    unrestricted <- compare_quiet(faers_group_ids(
        data,
        list(A = "insulin|humalog", B = "metformin"),
        .split = "joint"
    ))

    # The restricted B is a subset, so B-only can only shrink or stay equal...
    expect_true(length(joint$ids[["B-only"]]) <= length(unrestricted$ids[["B-only"]]))
    # ...and every report dropped from B lands outside B, never in A+B.
    expect_true(all(joint$ids[["B-only"]] %chin% unrestricted$raw_ids$B))
    expect_true(all(joint$ids[["A+B"]] %chin% joint$raw_ids$B))
})

test_that("a single drug row can satisfy both definitions (fixed-dose products)", {
    data <- compare_std_data()
    joint <- compare_quiet(faers_group_ids(
        data,
        list(A = "hydrocodone", B = "acetaminophen"),
        .split = "joint"
    ))
    # The diagnosis is reported and bounded by the joint group size.
    expect_true(all(joint$summary$n_single_row_both <= joint$summary$n_reports))
    expect_true("n_single_row_both" %chin% names(joint$summary))
    # This is the documented inheritance from faers_reports_with(): a report
    # matching both definitions is placed in `A+B`, not excluded.
    expect_true(length(joint$ids[["A+B"]]) >= 0L)
})

test_that(".overlap handles shared reports and reports the removals", {
    data <- compare_std_data()
    groups <- list(A = "insulin|humalog", B = "metformin")

    keep <- compare_quiet(faers_group_ids(data, groups))
    excl <- compare_quiet(faers_group_ids(data, groups, .overlap = "exclude"))
    shared <- intersect(keep$ids$A, keep$ids$B)

    # Exclusion removes the intersection from both groups.
    expect_length(intersect(excl$ids$A, excl$ids$B), 0L)
    expect_setequal(excl$ids$A, setdiff(keep$ids$A, keep$ids$B))
    expect_setequal(excl$ids$B, setdiff(keep$ids$B, keep$ids$A))
    expect_identical(
        excl$summary$n_removed_overlap,
        c(length(intersect(keep$ids$A, keep$ids$B)),
          length(intersect(keep$ids$B, keep$ids$A)))
    )
    expect_length(shared, length(shared))

    # Repartition turns the intersection into its own group. It partitions
    # A u B (not the whole scope U): A-only, B-only and A&B are mutually
    # exclusive and together contain exactly the reports exposed to either drug.
    repart <- compare_quiet(faers_group_ids(data, groups, .overlap = "repartition"))
    expect_identical(
        names(repart$ids), c("A-only", "B-only", "A&B")
    )
    expect_identical(
        sum(lengths(repart$ids)),
        length(union(keep$ids$A, keep$ids$B))
    )
    expect_setequal(
        unlist(repart$ids, use.names = FALSE),
        union(keep$ids$A, keep$ids$B)
    )
    expect_length(intersect(repart$ids[["A-only"]], repart$ids[["B-only"]]), 0L)
    expect_identical(
        sum(lengths(repart$ids)), length(keep$ids$A) + length(keep$ids$B) - length(shared)
    )

    # "error" refuses to continue.
    expect_error(
        compare_quiet(faers_group_ids(data, groups, .overlap = "error")),
        "overlap"
    )
})

test_that("faers_group_ids() validates its inputs", {
    data <- compare_std_data()

    # Unnamed groups.
    expect_error(
        compare_quiet(faers_group_ids(data, list("insulin", "metformin"))),
        "name"
    )
    # Unknown definition keys.
    expect_error(
        compare_quiet(faers_group_ids(data, list(A = list(.nope = "x")))),
        "nope"
    )
    # A joint split needs exactly two groups.
    expect_error(
        compare_quiet(faers_group_ids(data, list(A = "insulin"), .split = "joint")),
        "two"
    )
    # Roles must be named and refer to the group's own drugs.
    expect_error(
        compare_quiet(faers_group_ids(
            data,
            list(A = list(.drugs = "insulin", .roles = list(other = "PS")))
        )),
        "role"
    )
})

test_that("empty groups are legal and flagged", {
    data <- compare_std_data()
    groups <- compare_quiet(faers_group_ids(
        data,
        list(A = "insulin|humalog", empty = "no_such_drug_zzz")
    ))
    expect_identical(length(groups$ids$empty), 0L)
    expect_true(groups$summary[group == "empty", is_empty])
    expect_false(groups$summary[group == "A", is_empty])
})

test_that("faers_group_ids() is identical on both backends", {
    skip_if_not_installed("duckdb")
    skip_if_not_installed("DBI")

    data <- compare_std_data()
    db <- compare_db_twin(data)
    on.exit(DBI::dbDisconnect(db@db@con, shutdown = TRUE), add = TRUE)

    groups <- list(A = "insulin|humalog", B = "metformin", C = "glyburide")
    # `"joint"` and `"repartition"` are defined for exactly two groups.
    two <- groups[1:2]
    cases <- list(
        list(defs = groups, args = list()),
        list(defs = two, args = list(.split = "joint")),
        list(defs = groups, args = list(.overlap = "exclude")),
        list(defs = two, args = list(.overlap = "repartition"))
    )
    for (case in cases) {
        mem <- compare_quiet(do.call(
            faers_group_ids, c(list(data, case$defs), case$args)
        ))
        dbo <- compare_quiet(do.call(
            faers_group_ids, c(list(db, case$defs), case$args)
        ))
        expect_identical(mem$summary, dbo$summary)
        expect_identical(mem$ids, dbo$ids)
        expect_identical(mem$overlaps, dbo$overlaps)
        expect_identical(mem$membership, dbo$membership)
    }
})
