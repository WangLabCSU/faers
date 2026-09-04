# Integration tests for the optional DuckDB backend (`database = "duckdb"`).
#
# These tests are skipped when duckdb is not installed, so the default
# (`database = "memory"`) path and the rest of the test suite never depend on
# it.  They assert byte-parity between memory and duckdb for every
# *deterministic* output of the analysis pipeline.  The phv_signal MCMC
# columns are inherently stochastic (they differ run-to-run even on the memory
# path) and are therefore NOT compared byte-for-byte here.

has_duckdb <- requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("DBI", quietly = TRUE)

# Build the sample FAERSascii with all standardization state intact.
sample_object <- function() {
    rds <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
    sl <- list()
    for (sn in methods::slotNames("FAERSascii")) {
        sl[[sn]] <- tryCatch(methods::slot(rds, sn), error = function(e) NULL)
    }
    do.call(methods::new, c(list(Class = "FAERSascii"), sl))
}

# Build a db-backed twin of a memory object by ingesting each field table into
# a fresh in-memory duckdb.
db_twin <- function(mem_obj) {
    ns <- asNamespace("faers")
    con <- ns$db_connect(":memory:")
    tables <- ns$db_field_tables()
    for (f in names(mem_obj@data)) ns$db_ingest(con, f, mem_obj@data[[f]])
    ns$db_ingest_meddra(con, mem_obj)
    proxies <- stats::setNames(
        lapply(names(mem_obj@data), function(f) methods::new("FAERSdbTbl", con = con, table = tables[[f]])),
        names(mem_obj@data)
    )
    methods::new("FAERSascii",
        data = proxies,
        deletedCases = methods::slot(mem_obj, "deletedCases"),
        year = mem_obj@year, quarter = mem_obj@quarter,
        standardization = mem_obj@standardization,
        deduplication = mem_obj@deduplication, format = mem_obj@format,
        meddra = mem_obj@meddra,
        db = methods::new("FAERSdb", con = con, path = ":memory:",
            version = as.character(packageVersion("duckdb")), tables = tables)
    )
}

testthat::test_that("duckdb mode preserves parse pipeline", {
    testthat::skip_if_not(has_duckdb)
    dir <- system.file("extdata", package = "faers")
    mem <- faers(c(2004, 2017), c("q1", "q2"), "ascii",
        dir = dir, compress_dir = tempdir())
    db <- faers(c(2004, 2017), c("q1", "q2"), "ascii",
        dir = dir, compress_dir = tempdir(), database = "duckdb")

    testthat::expect_false(is.null(db@db))
    testthat::expect_s4_class(db@data$demo, "FAERSdbTbl")
    testthat::expect_true(data.table::is.data.table(faers_get(mem, "reac")))
    testthat::expect_true(data.table::is.data.table(faers_get(db, "reac")))
    testthat::expect_equal(faers_get(mem, "reac"), faers_get(db, "reac"))

    mat <- faers:::materialize(db)
    testthat::expect_null(mat@db)
    for (f in names(mem@data)) {
        testthat::expect_equal(mem@data[[f]], mat@data[[f]])
    }
})

testthat::test_that("faers_primaryid matches in db mode", {
    testthat::skip_if_not(has_duckdb)
    mem <- sample_object()
    db <- db_twin(mem)
    testthat::expect_identical(faers_primaryid(mem), faers_primaryid(db))
})

testthat::test_that("faers_counts matches in db mode (single/multi event, na.rm, fn)", {
    testthat::skip_if_not(has_duckdb)
    mem <- sample_object()
    db <- db_twin(mem)

    testthat::expect_identical(
        faers_counts(mem, .events = "soc_name", .field = "reac"),
        faers_counts(db, .events = "soc_name", .field = "reac")
    )
    testthat::expect_identical(
        faers_counts(mem, .events = c("soc_name", "soc_code"), .field = "reac"),
        faers_counts(db, .events = c("soc_name", "soc_code"), .field = "reac")
    )
    testthat::expect_identical(
        faers_counts(mem, .events = "soc_name", .field = "reac", .na.rm = TRUE),
        faers_counts(db, .events = "soc_name", .field = "reac", .na.rm = TRUE)
    )
    testthat::expect_identical(
        faers_counts(mem, .events = "soc_name", .field = "reac",
            .fn = function(d) d[!is.na(d$soc_name)]),
        faers_counts(db, .events = "soc_name", .field = "reac",
            .fn = function(d) d[!is.na(d$soc_name)])
    )
})

testthat::test_that("faers_phv_table contingency table matches in db mode", {
    testthat::skip_if_not(has_duckdb)
    mem <- sample_object()
    db <- db_twin(mem)
    testthat::expect_identical(
        suppressWarnings(faers_phv_table(mem, .full = mem, .events = "soc_name", .field = "reac")),
        suppressWarnings(faers_phv_table(db, .full = db, .events = "soc_name", .field = "reac"))
    )
})

testthat::test_that("faers_data and [[ materialize db objects identically", {
    testthat::skip_if_not(has_duckdb)
    mem <- sample_object()
    db <- db_twin(mem)
    testthat::expect_equal(faers_data(mem), faers_data(db))
    testthat::expect_equal(mem[["reac"]], db[["reac"]])
    testthat::expect_equal(mem$reac, db$reac)
})

testthat::test_that("faers_filter keeps only matching primaryids in db mode (regression)", {
    # Regression: faers_keep_db had the DELETE operator inverted, so db mode
    # kept everything EXCEPT the filtered subset (200 -> 185 instead of 15).
    testthat::skip_if_not(has_duckdb)
    mem <- sample_object()
    db <- db_twin(mem)
    fn <- function(d) unique(d$primaryid[
        grepl("METFORMIN|ASPIRIN|HUMALOG|AVANDIA", d$drugname)])
    flt_mem <- faers_filter(mem, .field = "drug", .fn = fn)
    flt_db  <- faers_filter(db,  .field = "drug", .fn = fn)
    testthat::expect_identical(faers_primaryid(flt_mem), faers_primaryid(flt_db))
    # both should keep only the 15 matching reports, not 200 - 15 = 185
    testthat::expect_equal(nrow(flt_mem@data$demo),
        length(unique(fn(faers_get(mem, "drug")))))
    testthat::expect_equal(faers_primaryid(flt_db), faers_primaryid(flt_mem))
})

