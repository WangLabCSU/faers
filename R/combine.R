#' Combine FAERS objects from different Quarterly files.
#'
#' Packed all `FAERSascii` or `FAERSxml` objects into a single `FAERSascii` or
#' `FAERSxml` object. It is important to note that all data passed to these
#' functions via the `...` argument must belong to the different [FAERS]
#' objects, indicating that they have the different period data (as defined by
#' [faers_period]).
#'
#' @param ... Multiple `FAERSxml` or `FAERSascii` objects or a list containing
#' `FAERSxml` or `FAERSascii` objects. Objects can be standardized by
#' [faers_standardize] but cannot be de-duplicated by [faers_dedup]. If we
#' combine deduplicated objects from different quarterly data files, duplicate
#' reports will be introduced again.
#' @return A [FAERSxml] or [FAERSascii] object.
#' @examples
#' # the files included in the package are sampled
#' data1 <- faers_parse(
#'     system.file("extdata", "aers_ascii_2004q1.zip", package = "faers"),
#'     compress_dir = tempdir()
#' )
#' data2 <- faers_parse(
#'     system.file("extdata", "faers_ascii_2017q2.zip", package = "faers"),
#'     compress_dir = tempdir()
#' )
#' faers_combine(data1, data2)
#' @export
faers_combine <- function(...) {
    x <- wrap_faers(...)
    l <- length(x)
    if (l == 1L) {
        cli::cli_alert(
            "Nothing to do since only one {.cls FAERS} data provided"
        )
        return(x[[1L]])
    }
    cli::cli_alert("Combining all {l} {.cls FAERS} Data{?s}")
    combine_faers(x)
}

combine_faers <- function(x, call = rlang::caller_env()) {
    type <- combine_faers_list_type(x, call = call)
    period <- combine_faers_period(x, call = call)
    is_dedup <- combine_faers_deduplication(x, call = call)
    is_stand <- combine_faers_standardization(x, call = call)
    meddra <- combine_faers_meddra(x, standardization = is_stand, call = call)
    is_db <- vapply(x, function(obj) !is.null(obj@db), logical(1L))
    if (type == "ascii" && all(is_db)) {
        elts <- combine_faers_db(x)
        return(methods::new("FAERSascii",
            data = elts$data,
            deletedCases = unique(
                unlist(lapply(x, faers_deleted_cases), use.names = FALSE)
            ),
            year = period$year, quarter = period$quarter,
            standardization = is_stand, deduplication = is_dedup,
            meddra = meddra,
            db = elts$db
        ))
    }
    if (type == "ascii" && any(is_db)) {
        cli::cli_abort("Cannot mix {.arg database = \"memory\"} and {.arg database = \"duckdb\"} objects when combining.")
    }
    data <- combine_faers_data(x, type)
    switch(type,
        ascii = methods::new("FAERSascii",
            data = data,
            deletedCases = unique(
                unlist(lapply(x, faers_deleted_cases), use.names = FALSE)
            ),
            year = period$year, quarter = period$quarter,
            standardization = is_stand, deduplication = is_dedup,
            meddra = meddra
        ),
        xml = methods::new("FAERSxml",
            data = data,
            year = period$year, quarter = period$quarter,
            standardization = is_stand, deduplication = is_dedup,
            meddra = meddra
        )
    )
}

combine_faers_list_type <- function(lst, call = rlang::caller_env()) {
    type <- NULL
    for (allowed_type in c("ascii", "xml")) {
        if (all(is_matched_faers(lst, allowed_type))) {
            type <- allowed_type
        }
    }
    if (is.null(type)) {
        cli::cli_abort(
            "All elements in {.arg ...} must be of the same class, {.cls FAERSascii} or {.cls FAERSxml}",
            call = call
        )
    }
    type
}

combine_faers_period <- function(lst, call = rlang::caller_env()) {
    out <- data.table::rbindlist(lapply(lst, faers_period))
    if (anyDuplicated(out)) {
        cli::cli_abort(c(
            "duplicated FAERS quarterly datas cannot be combined",
            i = "You can check {.fn faers_period}"
        ), call = call)
    }
    out
}

combine_faers_deduplication <- function(lst, call = rlang::caller_env()) {
    dedup_vec <- vapply(lst, function(object) object@deduplication, logical(1L))
    if (all(dedup_vec)) {
        # If we combine deduplicated objects from different quarterly data
        # files, duplicate reports will be introduced. Therefore, we should only
        # combine objects before performing deduplication.
        cli::cli_abort(
            "De-duplicated data must not be combined, you should always do de-duplication as a whole",
            call = call
        )
        TRUE
    } else if (!any(dedup_vec)) {
        FALSE
    } else {
        cli::cli_abort(
            "All elements in {.arg ...} must be fully undeduplicated.",
            call = call
        )
    }
}

combine_faers_standardization <- function(lst, call = rlang::caller_env()) {
    standardization_vec <- vapply(
        lst, function(object) object@standardization, logical(1L)
    )
    if (all(standardization_vec)) {
        TRUE
    } else if (!any(standardization_vec)) {
        FALSE
    } else {
        cli::cli_abort(
            "All elements in {.arg ...} must be either fully standardized or not at all.",
            call = call
        )
    }
}

combine_faers_data <- function(x, type) {
    is_db <- vapply(x, function(obj) !is.null(obj@db), logical(1L))
    if (type == "ascii" && all(is_db)) {
        return(combine_faers_db(x))
    }
    if (type == "ascii" && any(is_db)) {
        cli::cli_abort("Cannot mix {.arg database = \"memory\"} and {.arg database = \"duckdb\"} objects when combining.")
    }
    switch(type,
        ascii = combine_faers_ascii_data(x),
        xml = combine_faers_xml_data(x)
    )
}

combine_faers_ascii_data <- function(x) {
    data_list <- lapply(FAERS_ASCII_FILE_FIELDS, function(field) {
        data.table::rbindlist(
            lapply(x, function(obj) obj@data[[field]]),
            fill = TRUE, use.names = TRUE
        )
    })
    data.table::setattr(data_list, "names", FAERS_ASCII_FILE_FIELDS)
    data_list
}

combine_faers_xml_data <- function(x) {
    data.table::rbindlist(
        lapply(x, function(obj) obj@data),
        fill = TRUE, use.names = TRUE
    )
}

combine_faers_meddra <- function(lst, standardization, 
                                call = rlang::caller_env()) {
    if (!standardization) {
        return(NULL)
    }
    ref <- lst[[1L]]@meddra
    is_right <- vapply(
        lst, function(object) identical(object@meddra, ref), logical(1L)
    )
    if (all(is_right)) {
        return(ref)
    }
    cli::cli_abort(
        "All elements in {.arg ...} must use the same {.cls MedDRA} data",
        call = call
    )
}

is_matched_faers <- function(x, type) {
    vapply(x, methods::is, logical(1L), class2 = paste0("FAERS", type))
}

wrap_faers <- function(..., call = rlang::caller_env()) {
    if (!...length()) {
        cli::cli_abort("empty {.arg ...}", call = call)
    }
    if (is.list(..1)) {
        if (...length() == 1L) {
            ..1
        } else {
            cli::cli_abort(
                "You can only supply a single list or multiple {.cls FAERS} objects into {.arg ...}",
                call = call
            )
        }
    } else {
        list(...)
    }
}
