#' Counting the number of unique case for each event
#'
#' @param .object A [FAERSascii] object.
#' @return A [data.table][data.table::data.table] object.
#' @examples
#' # you must change `dir`, as the files included in the package are sampled
#' data <- faers(c(2004, 2017), c("q1", "q2"),
#'     dir = system.file("extdata", package = "faers"),
#'     compress_dir = tempdir()
#' )
#' \dontrun{
#' # you must standardize and deduplication before disproportionality analysis
#' # you should replace `meddra_path` with yours
#' data <- faers_standardize(data, meddra_path)
#' data <- faers_dedup(data)
#' faers_counts(data)
#' }
#' std_data <- readRDS(system.file("extdata", "standardized_data.rds",
#'     package = "faers"
#' ))
#' faers_counts(std_data)
#' @export
#' @name faers_counts
methods::setGeneric("faers_counts", function(.object, ...) {
    standardGeneric("faers_counts")
})

#' @param ... Other arguments passed to specific methods, for `FAERSascii`
#' method, other arguments passed to `.fn()`.
#' @param .field A string indicates the interested FAERS fields to
#' use. Only values "demo", "drug", "indi", "ther", "reac", "rpsr", and "outc"
#' can be used.
#' @param .events A character specify the events column(s) in the `.field` data
#' to count the unique `primaryid`. If multiple columns were selected, the
#' combination for all columns will define the interested events.
#' @param .fn A function or formula defined the preprocessing function before
#' creating contingency table, with the `.field` data as the input and return a
#' [data.table][data.table::data.table].
#'
#' Note: When using the `set*` or `:=` function from `data.table` with the
#' "demo", "drug", "ther", "rpsr", and "outc" data, exercise caution as these
#' functions directly modify the internal data. In such cases, it is advisable
#' to use the [copy][data.table::copy] function first.
#'
#'   If a **function**, it is used as is.
#'
#'   If a **formula**, e.g. `~ .x + 2`, it is converted to a function with up to
#'   two arguments: `.x` (single argument) or `.x` and `.y` (two arguments). The
#'   `.` placeholder can be used instead of `.x`.  This allows you to create
#'   very compact anonymous functions (lambdas) with up to two inputs.
#'
#'   If a **string**, the function is looked up in `globalenv()`.
#'
#' @param .na.rm A bool, whether `NA` value in `.events` column(s) should be
#' removed.
#' @rdname faers_counts
#' @export
methods::setMethod(
    "faers_counts", c(.object = "FAERSascii"),
    function(.object, .events = "soc_name", .fn = NULL, ..., .field = "reac", .na.rm = FALSE) {
        assert_bool(.na.rm)
        if (!.object@standardization) {
            cli::cli_abort("{.arg .object} must be standardized using {.fn faers_standardize}")
        }
        if (!is.null(.object@db)) {
            return(counts_db(.object, .events = .events, .fn = .fn, ...,
                .field = .field, .na.rm = .na.rm
            ))
        }
        data <- faers_get(.object, field = .field)
        if (!is.null(.fn)) {
            data <- rlang::as_function(.fn)(data, ...)
            if (!data.table::is.data.table(data)) {
                cli::cli_abort("{.fn .fn} must return a {.cls data.table}")
            }
        }
        groups <- c("primaryid", .events)
        data <- unique(data, by = groups, cols = character())
        if (.na.rm) {
            keep <- !Reduce(`|`, lapply(data[, .SD, .SDcols = .events], is.na))
            data <- data[keep]
        }
        eval(substitute(
            data[, list(N = .N), keyby = .events],
            list(.events = .events)
        ))
    }
)
