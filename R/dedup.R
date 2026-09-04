#' Tidy up FAERS Quarterly Data with duplicate records removed
#' @param object A [FAERSascii] object.
#' @param ... Other arguments passed to specific methods.
#' @return A [FAERSascii] object.
#' @seealso [faers_standardize]
#' @examples
#' # you must change `dir`, as the files included in the package are sampled
#' data <- faers(c(2004, 2017), c("q1", "q2"),
#'     dir = system.file("extdata", package = "faers"),
#'     compress_dir = tempdir()
#' )
#' \dontrun{
#' # we must standardize firstly
#' # you should replace `meddra_path` with yours
#' data <- faers_standardize(data, meddra_path)
#' faers_dedup(data)
#' }
#' @export
#' @name faers_dedup
methods::setGeneric("faers_dedup", function(object, ...) {
    standardGeneric("faers_dedup")
})

#' @param remove_deleted_cases If `TRUE`, will remove all
#' [deletedCases][FAERS-class] from the final result.
#' @export
#' @method faers_dedup FAERSascii
#' @rdname faers_dedup
methods::setMethod("faers_dedup", "FAERSascii", function(object, remove_deleted_cases = TRUE) {
    if (object@deduplication) {
        cli::cli_abort("You must not run {.fn faers_dedup} twice")
    }
    if (!object@standardization) {
        cli::cli_abort("{.cls FAERS} object must be standardized using {.fn faers_standardize} firstly")
    }
    if (!is.null(object@db)) {
        return(dedup_db(object, remove_deleted_cases = remove_deleted_cases))
    }
    if (isTRUE(remove_deleted_cases)) {
        deleted_cases <- faers_deleted_cases(object)
        if (!length(deleted_cases)) {
            deleted_cases <- NULL
        }
    } else {
        deleted_cases <- NULL
    }
    deduplicated_data <- do.call(
        dedup_faers_ascii,
        list(
            data = faers_mget(
                object,
                c("demo", "drug", "indi", "ther", "reac")
            ),
            deleted_cases = deleted_cases
        )
    )
    object@data <- lapply(object@data, function(x) {
        x[deduplicated_data,
            on = c("year", "quarter", "primaryid"),
            nomatch = NULL
        ]
    })
    object@deduplication <- TRUE
    object
})

#' @export
#' @method faers_dedup FAERSxml
#' @rdname faers_dedup
methods::setMethod("faers_dedup", "FAERSxml", function(object) {
    cli::cli_abort("Don't implement currently")
})

#' @method faers_dedup ANY
#' @rdname faers_dedup
methods::setMethod("faers_dedup", "ANY", function(object) {
    cli::cli_abort("Only {.cls FAERSascii} can work")
})

# define_common_by <- function(x) {
#     out <- NULL
#     nms <- names(x)
#     out <- c(out, choose_one(nms, c("primaryid", "isr")))
#     out <- c(out, choose_one(nms, c("caseid", "case")))
#     out
# }
# choose_one <- function(all, choices) {
#     for (choice in choices) {
#         if (any(choice == all)) {
#             return(choice)
#         }
#     }
#     cli::cli_abort("One of {.val {choices}} must exist")
# }

# Combine data from various source tables
# demo
# drug
# indi
# ther
# reac
# Perform string aggregation operations

dedup_faers_ascii <- function(data, deleted_cases = NULL) {
    # As recommended by the FDA, a deduplication step was performed to retain
    # the most recent report for each case with the same case identifier
    # nolint start
    cli::cli_alert("deduplication from the same source by retain the most recent report")

    # While the previous Legacy AERS (LAERS) database was Individual Safety
    # Report (ISR) based, the new FAERS database is Case/Version based. In
    # LAERS, a Case consisted of one or more ISRs (Initial and Follow-up
    # reports), and each ISR number represented a separate version of a case. A
    # case could contain multiple ISRs, and the latest ISR represented the most
    # current information about a particular case. (For example, Follow-up 1
    # would have the most up-to-date information about a case containing both an
    # Initial ISR and a Follow-up 1 ISR).
    # dt <- data.table::as.data.table(mtcars)
    # bench::mark(
    #   head = dt[order(mpg, gear), head(.SD, 1L), cyl],
    #   SD = dt[order(mpg, gear), .SD[1L], cyl],
    #   I = dt[dt[, .I[1L], cyl]$V1], # should be order firstly
    #   unique = unique(dt[order(mpg, gear)], by = "cyl"),
    #   unique2 = dt[order(mpg, gear), unique(.SD, by = "cyl")],
    #   check = FALSE, iterations = 100L
    # )
    #  expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc
    #   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int> <dbl>
    # 1 head        645.4µs  763.7µs      903.    66.6KB     9.12    99     1
    # 2 SD          648.8µs  721.7µs     1344.    66.6KB    27.4     98
    # 3 I             408µs    476µs     2012.      50KB     0      100
    # 4 unique       64.1µs   76.2µs    12337.    36.1KB   125.      99
    # 5 unique2     165.8µs  177.6µs     5363.    52.6KB     0      100

    # There are also duplicate reports where the same report was submitted by a
    # consumer and by the sponsor. we have found some duplicated `primaryid`
    # with different caseid in `demo`, so we remove this by keeping the latest
    # one, we remove `deleted_cases` firstly if it exist

    # 1. fda_dt: Date FDA received Case. In subsequent versions of a case, the
    #    latest manufacturer received date will be provided (YYYYMMDD format)
    # 2. i_f_code: I Initial; F Follow-up. latest should be "F". "F" < "I"
    # 3. event_dt: Date the adverse event occurred or began
    out <- data$demo
    if (!is.null(deleted_cases)) {
        out <- out[!deleted_cases, on = "caseid"]
    }
    data.table::setorderv(
        out,
        c("year", "quarter", "fda_dt", "i_f_code", "event_dt"),
        order = c(-1L, -1L, -1L, 1L, -1L),
        na.last = TRUE
    )
    out <- unique(out, by = "primaryid", cols = c(
        "year", "quarter", "caseid", "caseversion",
        "fda_dt", "i_f_code", "age_in_years", "sex",
        "country_code", "event_dt"
    ))

    # then we keep the `latest` informations for the patients
    # Such as caseid "11232882" in 2017q2 2019q2, 2019q3
    data.table::setorderv(out,
        cols = c(
            "year", "quarter", "caseversion", "fda_dt", "i_f_code", "event_dt"
        ),
        order = c(-1L, -1L, -1L, -1L, 1L, -1L), na.last = TRUE
    )
    out <- unique(out, by = "caseid")
    # test code
    # Following codes are based on
    # https://stackoverflow.com/questions/69366291/copy-only-one-variable-from-one-r-data-table-to-another-after-matching-on-a-vari
    # which should be memory efficient
    # a <- data.table::data.table(
    #     id = rep(letters[1:2], 5),
    #     var1 = rnorm(10)
    # )
    # b <- data.table::data.table(
    #     id = letters,
    #     var2 = rnorm(length(letters))
    # )
    # b[a[, paste0(var1, collapse = "/"), by = "id"],
    #     var1 := V1,
    #     on = "id"
    # ]
    # b
    # collapse all used drugs, indi, ther states, use it as a whole to identify
    # same cases.
    # match drug, indi, and ther data.
    common_keys <- c("year", "quarter", "primaryid")
    cli::cli_alert("merging `drug`, `indi`, `ther`, and `reac` data")
    out[
        data$drug[order(drug_seq),
            list(aligned_drugs = paste0(drugname, collapse = "/")),
            keyby = common_keys
        ],
        aligned_drugs := i.aligned_drugs,
        on = common_keys
    ]

    # should we remove unknown indications or just translate unknown indications
    # into NA ?
    # meddra_code: indi_pt
    # pt: 10070592 Product used for unknown indication
    # llt: 10057097 Drug use for unknown indication
    out[
        data$indi[order(indi_drug_seq, meddra_code),
            list(aligned_indi = paste0(meddra_code, collapse = "/")),
            keyby = common_keys
        ],
        aligned_indi := i.aligned_indi,
        on = common_keys
    ]
    out[
        data$ther[order(dsg_drug_seq, start_dt),
            list(aligned_start_dt = paste0(start_dt, collapse = "/")),
            keyby = common_keys
        ],
        aligned_start_dt := i.aligned_start_dt,
        on = common_keys
    ]
    # meddra_code: pt
    out[
        data$reac[order(meddra_code),
            list(aligned_reac = paste0(meddra_code, collapse = "/")),
            keyby = common_keys
        ],
        aligned_reac := i.aligned_reac,
        on = common_keys
    ]

    # consider two cases to be the same if they had a complete match of the
    # eight criteria which are sex, age, reporting country, event date, start
    # date, drug indications, drugs administered, and adverse reactions.  Two
    # records were also considered duplicated if they mismatch in only one of
    # the sex, age, reporting country, event date, start date, or drug
    # indications fields, but not the drug or adverse event fields.

    # Notes: always remember NA value in data.table, will be regarded as equal.
    # but they won't really be the same, so we should convert NA value into
    # other differentiated values, like ..__na_null__..1, ..__na_null__..2, and
    # so on.
    # round age_in_years to prevent minimal differences in age
    cli::cli_alert("deduplication from multiple sources by matching sex, age, reporting country, event date, start date, drug indications, drugs administered, and adverse reactions")
    out[, age_in_years_round := round(age_in_years, 2L)]
    can_be_ignored_columns <- c(
        "event_dt", "sex", "age_in_years_round", "country_code",
        "aligned_start_dt", "aligned_indi"
    )
    must_matched_columns <- c("aligned_drugs", "aligned_reac")
    all_columns <- c(must_matched_columns, can_be_ignored_columns)
    out[, (all_columns) := lapply(.SD, function(x) {
        idx <- is.na(x) | x == "NA"
        if (any(idx)) {
            x[idx] <- paste0("..__na_null__..", seq_len(sum(idx)))
        }
        x
    }), .SDcols = all_columns]
    # bench results:
    # dt <- data.table::as.data.table(mtcars)
    # dt2 <- data.table::as.data.table(mtcars)
    # dt3 <- data.table::as.data.table(mtcars)
    # keyed_dt <- data.table::as.data.table(mtcars)
    # data.table::setkeyv(keyed_dt, c("vs", "am"))
    # keyed_dt2 <- data.table::as.data.table(mtcars)
    # data.table::setkeyv(keyed_dt2, c("vs", "am"))
    # bench::mark(
    #   unique = unique(dt[order(-mpg, -gear)], by = "cyl"),
    #   unique2 = {
    #     data.table::setorderv(dt2, c("mpg", "gear"), order = -1L)
    #     unique(dt2, by = "cyl")
    #   },
    #   unique3 = {
    #     data.table::setkeyv(dt3, c("mpg", "gear"))
    #     unique(dt3, fromLast = TRUE, by = "cyl")
    #   },
    #   unique4 = unique(keyed_dt[order(-mpg, -gear)], by = "cyl"),
    #   unique5 = {
    #     data.table::setorderv(keyed_dt2, c("mpg", "gear"), order = -1L)
    #     unique(keyed_dt2, by = "cyl")
    #   },
    #   unique5 = dt[order(-mpg, -gear), unique(.SD, by = "cyl")],
    #   iterations = 100L, check = data.table::fsetequal
    # )
    # A tibble: 6 × 13
    #   expression      min  median `itr/sec` mem_alloc `gc/sec` n_itr
    #   <bch:expr> <bch:tm> <bch:t>     <dbl> <bch:byt>    <dbl> <int>
    # 1 unique       74.1µs    81µs    11974.    36.1KB      0     100
    # 2 unique2      31.1µs  33.7µs    28466.    16.3KB      0     100
    # 3 unique3      21.1µs  23.9µs    39440.    16.3KB      0     100
    # 4 unique4      71.7µs  76.9µs    11943.    36.1KB    121.     99
    # 5 unique5      29.3µs  30.8µs    30956.    16.6KB      0     100
    # 6 unique5     191.8µs 203.7µs     4704.    52.6KB     47.5    99

    for (i in seq_along(can_be_ignored_columns)) {
        # NAs are always first for `setkeyv`
        data.table::setkeyv(out, cols = common_keys)
        out <- unique(out,
            fromLast = TRUE,
            by = c(must_matched_columns, can_be_ignored_columns[-i])
        )
    }
    # out[, age_in_years_round := NULL]
    # nolint end
    # out[, (all_columns) := lapply(.SD, function(x) {
    #     x[startsWith(x, "..__na_null__..")] <- NA
    #     x
    # }), .SDcols = all_columns]
    # "caseversion", "fda_dt", "i_f_code",
    # "event_dt", "sex", "age_in_years", "country_code",
    # "aligned_drugs", "aligned_reac", "aligned_start_dt", "aligned_indi"
    out[, c("year", "quarter", "primaryid")]
}

utils::globalVariables(c(
    "drug_seq", "drugname", "indi_meddra_code", "start_dt",
    "indi_drug_seq", "dsg_drug_seq", "primaryid", "caseversion", "fda_dt",
    "i_f_code", "event_dt", "year", "caseid", "age_in_years_round",
    "meddra_code",
    paste0(c("", "i."), rep(
        c(
            "aligned_drugs", "aligned_indi", "aligned_start_dt",
            "aligned_reac"
        ),
        each = 2L
    ))
))
