testthat::test_that("Parsing external FAERS metadata works well", {
    testthat::skip_on_cran()
    testthat::skip_if_offline()
    testthat::expect_message(gds <- faers_meta(force = TRUE), "Reading html")
    testthat::expect_s3_class(gds, "data.table")
    testthat::expect_true(
        all(names(gds) == c("year", "quarter", "period", "ascii_urls", "ascii_file_size", "xml_urls", "xml_file_size"))
    )
    testthat::expect_true(all(grepl("^20\\d{2}$", gds$year)))
    testthat::expect_true(all(gds$quarter %in% paste0("q", 1:4)))
    testthat::expect_true(all(grepl("^(October|July|April|January)", gds$period)))
})

testthat::test_that("Parsing internal FAERS metadata works well", {
    if (file.exists(faers_meta_cache_file())) {
        file.remove(faers_meta_cache_file())
    }
    faers_clearcache("metadata")
    testthat::expect_message(
        gds <- faers_meta(force = FALSE, internal = TRUE), "internal"
    )
    testthat::expect_s3_class(gds, "data.table")
    testthat::expect_true(
        all(names(gds) == c("year", "quarter", "period", "ascii_urls", "ascii_file_size", "xml_urls", "xml_file_size"))
    )
    testthat::expect_true(all(grepl("^20\\d{2}$", gds$year)))
    testthat::expect_true(all(gds$quarter %in% paste0("q", 1:4)))
    testthat::expect_true(all(grepl("^(October|July|April|January)", gds$period)))
})
