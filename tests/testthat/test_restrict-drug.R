# Tests for faers_reports_with (drug-exposure report set constructor).
#
# The bundled sample (200 standardized reports) supplies realistic drug rows;
# a hand-crafted synthetic object covers exact AND/OR/role/exclude semantics,
# fixed-dose same-row matching and regex edge cases. duckdb parity tests
# reuse the twin-building pattern of test_duckdb.R / test_scan.R and are
# skipped when duckdb is not installed.

has_duckdb <- requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("DBI", quietly = TRUE)

reports_sample_object <- function() {
    rds <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
    sl <- list()
    for (sn in methods::slotNames("FAERSascii")) {
        sl[[sn]] <- tryCatch(methods::slot(rds, sn), error = function(e) NULL)
    }
    obj <- do.call(methods::new, c(list(Class = "FAERSascii"), sl))
    obj@deduplication <- TRUE
    obj
}

# A tiny standardized object with known drug contents:
#   1 -> ASPIRIN(PS), TYLENOL(C)   2 -> "aspirin "(SS)  [trailing space]
#   3 -> BUFFERIN(C)               4 -> IBUPROFEN(PS), ASPIRIN(SS)
#   5 -> "aspirin and dipyridamole" (one fixed-dose combination row)
reports_synth_object <- function() {
    obj <- reports_sample_object()
    obj@data$demo <- data.table::data.table(
        year = 2004L, quarter = "q1", primaryid = as.character(1:5)
    )
    obj@data$drug <- data.table::data.table(
        year = 2004L, quarter = "q1",
        primaryid = c("1", "1", "2", "3", "4", "4", "5"),
        drug_seq = 1:7,
        role_cod = c("PS", "C", "SS", "C", "PS", "SS", "PS"),
        drugname = c(
            "ASPIRIN", "TYLENOL", "aspirin ", "BUFFERIN",
            "IBUPROFEN", "ASPIRIN", "ASPIRIN AND DIPYRIDAMOLE"
        ),
        prod_ai = c(
            "ASPIRIN", "ACETAMINOPHEN", "ASPIRIN", "BUFFERIN",
            "IBUPROFEN", "ASPIRIN", "ASPIRIN AND DIPYRIDAMOLE"
        )
    )
    obj
}

testthat::test_that("faers_reports_with validates inputs", {
    obj <- reports_sample_object()

    # .drugs type/emptiness
    expect_error(faers_reports_with(obj, NULL), "character vector")
    expect_error(faers_reports_with(obj, character()), "character vector")
    expect_error(faers_reports_with(obj, c("a", "")), "empty strings")

    # duplicated patterns
    expect_error(faers_reports_with(obj, c("a", "a")), "duplicated")

    # invalid regex aborts with the offending pattern named
    expect_error(faers_reports_with(obj, "warfarin("), "Invalid regular expression")

    # bad roles: unknown drug, unknown code, empty, unnamed, duplicated name
    expect_error(
        faers_reports_with(obj, "aspirin", .roles = list(nope = "PS")),
        "not listed"
    )
    expect_error(
        faers_reports_with(obj, "aspirin", .roles = list(aspirin = "XX")),
        "Invalid role code"
    )
    expect_error(
        faers_reports_with(obj, "aspirin", .roles = list(aspirin = character())),
        "non-empty"
    )
    expect_error(faers_reports_with(obj, "aspirin", .roles = list("PS")), "named list")

    # standardization guard
    raw <- reports_sample_object()
    raw@standardization <- FALSE
    expect_error(faers_reports_with(raw, "aspirin"), "standardized")
})

testthat::test_that("faers_reports_with single pattern counts and sorts ids", {
    obj <- reports_synth_object()
    ids <- suppressMessages(faers_reports_with(obj, "aspirin", .return = "ids"))
    # reports 1, 2, 4, 5 (case/trailing-space variants included)
    expect_setequal(ids, c("1", "2", "4", "5"))
    expect_identical(as.character(ids), as.character(sort(ids)))
    expect_identical(
        suppressMessages(faers_reports_with(obj, "aspirin", .return = "count")),
        4L
    )
    # case-insensitive parity
    expect_identical(
        suppressMessages(faers_reports_with(obj, "ASPIRIN", .return = "ids")),
        ids
    )
})

testthat::test_that("faers_reports_with AND requires every pattern, OR any", {
    obj <- reports_synth_object()
    n_or <- suppressMessages(faers_reports_with(
        obj, c("aspirin", "ibuprofen"), .logic = "or", .return = "count"
    ))
    n_and <- suppressMessages(faers_reports_with(
        obj, c("aspirin", "ibuprofen"), .return = "count"
    ))
    expect_identical(n_or, 4L) # aspirin reports 1,2,4,5 + ibuprofen 4
    expect_identical(n_and, 1L) # only report 4 has both
})

testthat::test_that("faers_reports_with AND allows one row to satisfy two patterns", {
    obj <- reports_synth_object()
    # report 5's single row "ASPIRIN AND DIPYRIDAMOLE" must satisfy both
    ids <- suppressMessages(faers_reports_with(
        obj, c("aspirin", "dipyridamole"), .return = "ids"
    ))
    expect_true("5" %in% ids)
})

testthat::test_that("faers_reports_with roles are strict and effective", {
    obj <- reports_synth_object()
    # aspirin unrestricted: 4 reports; PS only: reports 1 (ASPIRIN PS)
    n_any <- suppressMessages(faers_reports_with(obj, "aspirin", .return = "count"))
    n_ps <- suppressMessages(faers_reports_with(
        obj, "aspirin", .roles = list(aspirin = "PS"), .return = "count"
    ))
    expect_identical(n_ps, 2L) # report 1 (PS) + report 4 (ASPIRIN PS row)
    expect_true(n_ps < n_any)
    # unspecified drugs stay unrestricted
    n_mixed <- suppressMessages(faers_reports_with(
        obj, c("aspirin", "ibuprofen"),
        .roles = list(aspirin = "PS"), .return = "count"
    ))
    expect_identical(n_mixed, 0L) # no report has both: ibuprofen AND aspirin-PS
    # OR with role restriction: aspirin PS (1,4) ∪ ibuprofen unrestricted (4)
    # ∪ dipyridamole (5, unrestricted) = 3
    n_or <- suppressMessages(faers_reports_with(
        obj, c("aspirin", "ibuprofen"),
        .roles = list(aspirin = "PS"), .logic = "or", .return = "count"
    ))
    expect_identical(n_or, 3L)
    expect_true(n_or <= 4L) # never exceeds unrestricted OR
})

testthat::test_that("faers_reports_with excludes reports at the report level", {
    obj <- reports_synth_object()
    # excluding aspirin also removes reports whose aspirin row is only C/SS
    ids <- suppressMessages(faers_reports_with(
        obj, "ibuprofen", .exclude_drugs = "aspirin", .return = "ids"
    ))
    expect_setequal(ids, character(0)) # report 4 also has aspirin
    # role-agnostic exclusion: C-role aspirin row also triggers exclusion
    n_with <- suppressMessages(faers_reports_with(obj, "ibuprofen", .return = "count"))
    n_without <- suppressMessages(faers_reports_with(
        obj, "ibuprofen", .exclude_drugs = "aspirin", .return = "count"
    ))
    expect_identical(n_without, 0L)
    expect_identical(n_with, 1L)
})

testthat::test_that("faers_reports_with handles regex special characters identically", {
    obj <- reports_synth_object()
    # dot + quantifier + anchors: "^aspi. n$" can never match (no "ASPIRN"),
    # a proven zero-hit; the class form "asp[i]rin" must match like grepl
    expect_identical(
        suppressMessages(faers_reports_with(obj, "^aspi. n$", .return = "ids")),
        character(0)
    )
    # character class (grepl semantics is the frozen baseline):
    # "asp[i]rin" matches ASPIRIN rows like grepl would
    ids_class <- suppressMessages(faers_reports_with(
        obj, "asp[i]rin", .return = "ids"
    ))
    expect_setequal(ids_class, c("1", "2", "4", "5"))
    # a literal '-' pattern still matches names containing hyphens
    ids_dash <- suppressMessages(faers_reports_with(obj, "AND DIPYRIDAMOLE", .return = "ids"))
    expect_setequal(ids_dash, c("5"))
})

testthat::test_that("faers_reports_with returns empty sets without error", {
    obj <- reports_synth_object()
    expect_identical(
        suppressMessages(faers_reports_with(obj, "zzzz", .return = "count")), 0L
    )
    expect_identical(
        suppressMessages(faers_reports_with(obj, "zzzznope", .return = "ids")),
        character(0)
    )
    o0 <- suppressMessages(faers_reports_with(obj, "zzzznope"))
    expect_s4_class(o0, "FAERSascii")
    expect_identical(nrow(faers_get(o0, "demo")), 0L)
    # OR with one zero-hit pattern still returns the other's hits
    n <- suppressMessages(faers_reports_with(
        obj, c("aspirin", "zzzznope"), .logic = "or", .return = "count"
    ))
    expect_identical(n, 4L)
})

testthat::test_that("faers_reports_with returns object, ids and count forms", {
    obj <- reports_synth_object()
    ids <- suppressMessages(faers_reports_with(obj, "aspirin", .return = "ids"))
    expect_type(ids, "character")
    expect_true(!anyDuplicated(ids))

    count <- suppressMessages(faers_reports_with(obj, "aspirin", .return = "count"))
    expect_identical(count, as.integer(length(ids)))

    o <- suppressMessages(faers_reports_with(obj, "aspirin"))
    expect_s4_class(o, "FAERSascii")
    expect_setequal(as.character(faers_get(o, "demo")$primaryid), ids)
    # flags survive
    expect_true(o@standardization)
    expect_true(o@deduplication)
    # input untouched
    expect_identical(nrow(faers_get(obj, "demo")), 5L)
    expect_identical(nrow(faers_get(obj, "drug")), 7L)
})

testthat::test_that("faers_reports_with pick-column logic", {
    obj <- reports_synth_object()
    # prod_ai matches too (same rows here), but drugname wins by default
    expect_silent(suppressMessages(
        faers_reports_with(obj, "dipyridamole", .return = "ids")
    ))
    # explicit drugname column: "bufferin" matches BUFFERIN (ignore.case)
    n <- suppressMessages(faers_reports_with(
        obj, "bufferin", .drug_col = "drugname", .return = "count"
    ))
    expect_identical(n, 1L)
    # unknown column aborts
    expect_error(
        suppressMessages(faers_reports_with(obj, "aspirin", .drug_col = "nope")),
        "not found"
    )
})

testthat::test_that("faers_reports_with missing name columns aborts clearly", {
    obj <- reports_synth_object()
    obj@data$drug$drugname <- NULL
    obj@data$drug$prod_ai <- NULL
    expect_error(suppressMessages(faers_reports_with(obj, "aspirin")), "explicitly")
})

testthat::test_that("faers_reports_with duplicated drug rows do not double count", {
    obj <- reports_synth_object()
    drug <- data.table::copy(faers_get(obj, "drug"))
    extra <- data.table::copy(drug[1L])
    extra$drug_seq <- 99L
    obj@data$drug <- data.table::rbindlist(list(drug, extra), use.names = TRUE)
    n <- suppressMessages(faers_reports_with(obj, "aspirin", .return = "count"))
    expect_identical(n, 4L) # duplicate rows never double-count reports
})

# Build a duckdb twin of a memory object (test_scan.R / test_duckdb.R pattern).
# Caller owns the connection lifecycle: pass a fresh con and defer its shutdown.
reports_db_twin <- function(mem_obj, con) {
    ns <- asNamespace("faers")
    tables <- ns$db_field_tables()
    for (f in names(mem_obj@data)) ns$db_ingest(con, f, mem_obj@data[[f]])
    ns$db_ingest_meddra(con, mem_obj)
    proxies <- stats::setNames(lapply(names(mem_obj@data), function(f) {
        methods::new("FAERSdbTbl", con = con, table = tables[[f]])
    }), names(mem_obj@data))
    methods::new("FAERSascii",
        data = proxies, deletedCases = methods::slot(mem_obj, "deletedCases"),
        year = mem_obj@year, quarter = mem_obj@quarter,
        standardization = TRUE, deduplication = TRUE,
        format = mem_obj@format, meddra = mem_obj@meddra,
        db = methods::new("FAERSdb", con = con, path = ":memory:",
            version = as.character(packageVersion("duckdb")), tables = tables)
    )
}

# Every pair of rows belongs to a real report in demo. Reports 1 and 2
# satisfy AND across opposite character partitions; 3 and 4 have exclusions
# in the opposite partition. Report 7 checks Unicode folding, and report 8
# satisfies both inclusion patterns on one combination-product row.
reports_mixed_object <- function() {
    obj <- reports_synth_object()
    obj@data$demo <- data.table::data.table(
        year = 2004L, quarter = "q1", primaryid = as.character(1:8)
    )
    obj@data$drug <- data.table::data.table(
        primaryid = as.character(c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 1)),
        drugname = c(
            "ASPIRIN", "caf\u00e9 warfarin", "caf\u00e9 aspirin", "WARFARIN",
            "caf\u00e9 aspirin", "vitamin k", "ASPIRIN", "caf\u00e9 vitamin k",
            "ASPIRIN", "WARFARIN", "caf\u00e9 aspirin", "caf\u00e9 warfarin",
            "ASPIRIN", "\u212aetamine", "aspirin and warfarin", "ASPIRIN"
        ),
        role_cod = c("PS", "SS", "PS", "PS", "C", "I", "SS", "C",
                     "PS", "SS", "PS", "PS", "PS", "PS", "PS", "PS")
    )
    obj@data$drug[, prod_ai := drugname]
    obj
}

testthat::test_that("mixed-character report sets combine before AND/OR and exclusion", {
    testthat::skip_if_not(has_duckdb)
    mem <- reports_mixed_object()
    con <- asNamespace("faers")$db_connect(":memory:")
    withr::defer(DBI::dbDisconnect(con, shutdown = TRUE))
    db <- reports_db_twin(mem, con)
    cases <- list(
        list(args = list(.drugs = c("aspirin", "warfarin")), ids = c("1", "2", "5", "6", "8")),
        list(args = list(.drugs = c("aspirin", "warfarin"), .logic = "or"), ids = as.character(1:8)),
        list(args = list(.drugs = c("aspirin", "warfarin"), .roles = list(warfarin = "PS")), ids = c("2", "6", "8")),
        list(args = list(.drugs = c("aspirin", "warfarin"), .roles = list(aspirin = "PS", warfarin = "SS")), ids = c("1", "5")),
        list(args = list(.drugs = "aspirin", .exclude_drugs = "vitamin k"), ids = c("1", "2", "5", "6", "7", "8")),
        list(args = list(.drugs = c("aspirin", "warfarin"), .logic = "or", .exclude_drugs = "vitamin k"), ids = c("1", "2", "5", "6", "7", "8")),
        list(args = list(.drugs = "aspirin", .exclude_drugs = "vitamin [k]"), ids = c("1", "2", "5", "6", "7", "8")),
        list(args = list(.drugs = c("asp[i]rin", "warfarin")), ids = c("1", "2", "5", "6", "8")),
        list(args = list(.drugs = c("aspirin", "zzzzmissing")), ids = character()),
        list(args = list(.drugs = "aspirin", .exclude_drugs = "aspirin"), ids = character())
    )
    for (col in c("drugname", "prod_ai")) for (case in cases) {
        for (object in list(mem, db)) {
            args <- c(list(.object = object, .drug_col = col), case$args)
            actual <- suppressMessages(do.call(faers_reports_with, c(args, list(.return = "ids"))))
            expect_identical(actual, case$ids)
            expect_identical(
                suppressMessages(do.call(faers_reports_with, c(args, list(.return = "count")))),
                as.integer(length(case$ids))
            )
        }
    }
})

testthat::test_that("Unicode folding and hit diagnostics use complete per-pattern sets", {
    testthat::skip_if_not(has_duckdb)
    mem <- reports_mixed_object()
    con <- asNamespace("faers")$db_connect(":memory:")
    withr::defer(DBI::dbDisconnect(con, shutdown = TRUE))
    db <- reports_db_twin(mem, con)
    for (pattern in c("ketamine", "k", "s")) {
        expected <- sort(unique(mem@data$drug$primaryid[
            grepl(pattern, mem@data$drug$drugname, ignore.case = TRUE)
        ]))
        expect_identical(
            suppressMessages(faers_reports_with(db, pattern, .return = "ids")), expected
        )
        expected_exclusion <- sort(setdiff(as.character(1:8), expected))
        expect_identical(suppressMessages(faers_reports_with(
            db, "aspirin", .exclude_drugs = pattern, .return = "ids"
        )), expected_exclusion)
    }
    # Duplicate ASCII rows and cross-partition hits count once per report.
    expect_message(
        n <- faers_reports_with(db, c("aspirin", "warfarin"), .return = "count"),
        "aspirin: 8"
    )
    expect_identical(n, 5L)
    expect_message(
        faers_reports_with(db, "warfarin", .roles = list(warfarin = "PS"), .return = "count"),
        "warfarin: 3"
    )
})

testthat::test_that("mixed-character object returns preserve the source and previous results", {
    testthat::skip_if_not(has_duckdb)
    mem <- reports_mixed_object()
    con <- asNamespace("faers")$db_connect(":memory:")
    withr::defer(DBI::dbDisconnect(con, shutdown = TRUE))
    db <- reports_db_twin(mem, con)
    original <- lapply(names(mem@data), function(f) db[[f]])
    both <- suppressMessages(faers_reports_with(db, c("aspirin", "warfarin")))
    withr::defer({
        DBI::dbDisconnect(both@db@con, shutdown = TRUE)
        unlink(c(both@db@path, paste0(both@db@path, ".wal")))
    })
    no_k <- suppressMessages(faers_reports_with(db, "aspirin", .exclude_drugs = "vitamin k"))
    withr::defer({
        DBI::dbDisconnect(no_k@db@con, shutdown = TRUE)
        unlink(c(no_k@db@path, paste0(no_k@db@path, ".wal")))
    })
    expect_setequal(faers_primaryid(both), c("1", "2", "5", "6", "8"))
    expect_setequal(faers_primaryid(no_k), c("1", "2", "5", "6", "7", "8"))
    expect_identical(both@standardization, db@standardization)
    expect_identical(both@deduplication, db@deduplication)
    for (i in seq_along(mem@data)) {
        f <- names(mem@data)[[i]]
        expect_equal(db[[f]], original[[i]])
        expect_true(all(both[[f]]$primaryid %in% faers_primaryid(both)))
        expect_true(all(no_k[[f]]$primaryid %in% faers_primaryid(no_k)))
    }
})

testthat::test_that("literal report queries do not materialize the full drug table", {
    testthat::skip_if_not(has_duckdb)
    con <- asNamespace("faers")$db_connect(":memory:")
    withr::defer(DBI::dbDisconnect(con, shutdown = TRUE))
    db <- reports_db_twin(reports_mixed_object(), con)
    testthat::local_mocked_bindings(
        db_collect = function(...) stop("Full-table collection is not needed"),
        .package = "faers"
    )
    expect_identical(suppressMessages(faers_reports_with(
        db, c("aspirin", "warfarin"), .return = "count"
    )), 5L)
    expect_identical(suppressMessages(faers_reports_with(
        db, "aspirin", .exclude_drugs = "vitamin k", .return = "ids"
    )), c("1", "2", "5", "6", "7", "8"))
})

testthat::test_that("faers_reports_with parity on db twin (AND/OR/roles/exclude/regex)", {
    testthat::skip_if_not(has_duckdb)
    ns <- asNamespace("faers")
    mem <- reports_synth_object()
    con <- ns$db_connect(":memory:")
    withr::defer(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE))
    db_obj <- reports_db_twin(mem, con)

    # literal pattern: pushed down on db, plain grepl on memory
    expect_identical(
        suppressMessages(faers_reports_with(mem, "aspirin", .return = "ids")),
        suppressMessages(faers_reports_with(db_obj, "aspirin", .return = "ids"))
    )
    # OR + roles + exclude (pushed-down path)
    m <- suppressMessages(faers_reports_with(
        mem, c("aspirin", "ibuprofen"), .logic = "or",
        .roles = list(aspirin = "PS"), .exclude_drugs = "tylenol",
        .return = "ids"
    ))
    d <- suppressMessages(faers_reports_with(
        db_obj, c("aspirin", "ibuprofen"), .logic = "or",
        .roles = list(aspirin = "PS"), .exclude_drugs = "tylenol",
        .return = "ids"
    ))
    expect_identical(m, d)
    # regex pattern: falls back to R matching on db, identical to memory
    d_rx <- suppressMessages(faers_reports_with(db_obj, "asp[i]rin", .return = "ids"))
    expect_setequal(d_rx, c("1", "2", "4", "5"))
    # mixed pushdown/fallback also matches memory semantics
    d_mix <- suppressMessages(faers_reports_with(
        db_obj, c("aspirin", "^ibu"), .logic = "or", .return = "ids"
    ))
    expect_setequal(d_mix, c("1", "2", "4", "5"))
    m_mix <- suppressMessages(faers_reports_with(
        mem, c("aspirin", "^ibu"), .logic = "or", .return = "ids"
    ))
    expect_identical(d_mix, m_mix)
    expect_identical(
        suppressMessages(faers_reports_with(mem, "aspirin", .return = "count")),
        suppressMessages(faers_reports_with(db_obj, "aspirin", .return = "count"))
    )
})

testthat::test_that("faers_reports_with leaves the source db untouched across builds", {
    testthat::skip_if_not(has_duckdb)
    ns <- asNamespace("faers")
    mem <- reports_synth_object()
    con <- ns$db_connect(":memory:")
    withr::defer(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE))
    db_obj <- reports_db_twin(mem, con)
    n0 <- ns$db_nrow(con, "drug")
    # a sequence of constructors and object-mode returns on the SAME db
    o1 <- suppressMessages(faers_reports_with(db_obj, "aspirin"))
    o2 <- suppressMessages(faers_reports_with(db_obj, "ibuprofen"))
    i1 <- suppressMessages(faers_reports_with(db_obj, "aspirin", .return = "ids"))
    i2 <- suppressMessages(faers_reports_with(db_obj, "ibuprofen", .return = "ids"))
    expect_identical(ns$db_nrow(con, "drug"), n0) # untouched
    expect_identical(i1, c("1", "2", "4", "5"))
    expect_identical(i2, c("4"))
    # object-mode returns filter their own private copy, not the shared db
    expect_setequal(as.character(faers_get(o1, "demo")$primaryid), c("1", "2", "4", "5"))
    expect_setequal(as.character(faers_get(o2, "demo")$primaryid), c("4"))
    # and the original still sees the full data
    expect_identical(
        suppressMessages(faers_reports_with(db_obj, "aspirin", .return = "count")),
        4L
    )
})
