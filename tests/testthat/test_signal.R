testthat::test_that("`faers_phv_table()` works well", {
    testthat::skip_if_not(dir.exists("~/Data/MedDRA/MedDRA_26_1_English"))
    data <- faers(c(2004, 2017),
        c("q1", "q2"), "ascii",
        dir = internal_file("extdata"),
        compress_dir = tempdir()
    )
    data_std <- faers_standardize(data,
        "~/Data/MedDRA/MedDRA_26_1_English", # nolint
        add_smq = TRUE
    )
    testthat::expect_error(faers_phv_table(data))
    interested_ids <- sample(faers_primaryid(data_std), 100L)
    obj1 <- faers_keep(data_std, interested_ids)
    obj2 <- faers_keep(data_std, interested_ids, invert = TRUE)
    testthat::expect_no_error(phv_table <- faers_phv_table(
        faers_filter(data_std, .fn = function(x) {
            interested_ids
        }, .field = "demo"),
        .full = data_std
    ))
    phv_table1 <- faers_phv_table(obj1, .full = data_std)
    testthat::expect_true(data.table::fsetequal(phv_table, phv_table1))
    phv_table2 <- faers_phv_table(obj1, .object2 = obj2)
    testthat::expect_true(data.table::fsetequal(phv_table, phv_table2))
})

testthat::test_that("`faers_phv_signal()` works well", {
    testthat::skip_if_not(dir.exists("~/Data/MedDRA/MedDRA_26_1_English"))
    data <- faers(c(2004, 2017),
        c("q1", "q2"), "ascii",
        dir = internal_file("extdata"),
        compress_dir = tempdir()
    )
    data_std <- faers_standardize(data,
        "~/Data/MedDRA/MedDRA_26_1_English", # nolint
        add_smq = TRUE
    )
    interested_ids <- sample(faers_primaryid(data_std), 100L)
    obj1 <- faers_keep(data_std, interested_ids)
    obj2 <- faers_keep(data_std, interested_ids, invert = TRUE)
    testthat::expect_error(faers_phv_signal(obj1))
    testthat::expect_no_error(signal <- suppressWarnings(faers_phv_signal(
        faers_filter(data_std, .fn = function(x) {
            interested_ids
        }),
        .full = data_std,
        BPPARAM = BiocParallel::SerialParam(RNGseed = 1L)
    )))
    signal1 <- suppressWarnings(faers_phv_signal(obj1,
        .full = data_std,
        BPPARAM = BiocParallel::SerialParam(RNGseed = 1L)
    ))
    testthat::expect_true(data.table::fsetequal(signal, signal1))
    signal2 <- suppressWarnings(faers_phv_signal(obj1,
        .object2 = obj2,
        BPPARAM = BiocParallel::SerialParam(RNGseed = 1L)
    ))
    testthat::expect_true(data.table::fsetequal(signal, signal2))
})

# Regression: the vectorized methods that build their result from a list of
# per-row results (`chisq`, `fisher`, `bcpnn_mcmc`) used to fail on zero-length
# input, because `transpose()`/`rbind()` of an empty list yields a zero-column
# data.table and the following `setnames()` then had nothing to name. They now
# return their full schema, so a zero-length call stays bindable and keeps the
# same columns as a populated one.
testthat::test_that("phv_* methods keep their schema on zero-length input", {
    empty <- list(a = numeric(), b = numeric(), c = numeric(), d = numeric())
    methods <- c("ror", "prr", "chisq", "fisher", "bcpnn_norm", "obsexp_shrink", "ebgm")
    populated <- suppressWarnings(suppressMessages(do.call(
        phv_signal,
        c(list(a = c(5, 10, 3), b = c(2, 3, 0), c = c(1, 2, 4), d = c(100, 200, 50)),
            list(methods = methods))
    )))
    zero <- suppressWarnings(suppressMessages(do.call(
        phv_signal, c(empty, list(methods = methods))
    )))
    testthat::expect_identical(nrow(zero), 0L)
    testthat::expect_identical(names(zero), names(populated))
    testthat::expect_identical(
        vapply(zero, function(x) class(x)[[1L]], character(1L)),
        vapply(populated, function(x) class(x)[[1L]], character(1L))
    )

    # `chisq`/`fisher`/`bcpnn_mcmc` are the ones that regressed; assert them
    # individually so a future rewrite cannot silently drop a companion column
    # (e.g. `chisq` without `chisq_pvalue`).
    testthat::expect_identical(
        names(suppressWarnings(phv_chisq(numeric(), numeric(), numeric(), numeric()))),
        c("chisq", "pvalue")
    )
    testthat::expect_identical(
        names(suppressWarnings(phv_fisher(numeric(), numeric(), numeric(), numeric()))),
        c("odds_ratio", "ci_low", "ci_high", "pvalue")
    )
    testthat::expect_identical(
        names(suppressWarnings(
            phv_bcpnn_mcmc(numeric(), numeric(), numeric(), numeric(), n_mcmc = 100L)
        )),
        c("ic", "ci_low", "ci_high")
    )

    # MCMC-free methods are also callable with nothing to compute.
    testthat::expect_identical(nrow(phv_ror(numeric(), numeric(), numeric(), numeric())), 0L)
    testthat::expect_identical(nrow(phv_prr(numeric(), numeric(), numeric(), numeric())), 0L)
    testthat::expect_identical(nrow(phv_bcpnn_norm(numeric(), numeric(), numeric(), numeric())), 0L)

    # The companion issue: a non-empty call must still fail on malformed input,
    # so the guard cannot have swallowed the validation.
    testthat::expect_error(phv_chisq(-1, 1, 1, 1), "non-negative")
    testthat::expect_error(phv_fisher(1, 1, 1, 1:2), "same length")
    testthat::expect_error(phv_bcpnn_mcmc(1, 1, NA_real_, 1), "non-negative")
})
