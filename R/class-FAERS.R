#' FAERS class
#'
#' @description Provide a container for FAERS Quarterly Data file
#' @slot year An integer specifies the year information.
#' @slot quarter A string specifies the quarter information.
#' @slot data For `FAERSxml`, a [data.table][data.table::data.table]. For
#' `FAERSascii`, a list of [data.table][data.table::data.table].
#' @slot meddra A [MedDRA] or `NULL` representing the meddra data used for
#' standardization.
#' @slot format A string of "ascii" or "xml" indicates the file format used.
#' @slot deletedCases An atomic character, as of 2019 Quarter one there are new
#' files that lists deleted cases. [faers_dedup] will remove cases in this slot.
#' @slot standardization A bool, indicates whether standardization has been
#' performed.
#' @slot deduplication A bool, indicates whether deduplication has been
#' performed.
#' @details
#'  - `faers_data`: Extract the `data` slot.
#'  - `faers_year`: Extract the `year` slot.
#'  - `faers_quarter`: Extract the `quarter` slot.
#'  - `faers_period`: A [data.table][data.table::data.table] combine the `year`
#'    and `quarter` slot.
#'  - `faers_meddra`: Extract the `meddra` slot. If `object` have never been
#'    standardized, always return `NULL`.
#'  - `faers_deleted_cases`: Extract the `deletedCases` slot.
#' @return See details.
#' @examples
#' # ususaly we use faers() function to create a `FAERS` object
#' # you must change `dir`, as the file included in the package is sampled
#' data <- faers(2004, "q1",
#'     dir = system.file("extdata", package = "faers"),
#'     compress_dir = tempdir()
#' )
#' faers_data(data)
#' faers_year(data)
#' faers_quarter(data)
#' faers_period(data)
#' faers_meddra(data)
#' faers_deleted_cases(data)
#' @aliases FAERS
#' @name FAERS-class
NULL

#' @include meddra.R
#' @include class-FAERSdb.R
methods::setClassUnion("MedDRAOrNull", c("NULL", "MedDRA"))

#' @export
methods::setClass(
    "FAERS",
    slots = list(
        year = "integer",
        quarter = "character",
        data = "ANY",
        meddra = "MedDRAOrNull",
        deduplication = "logical",
        standardization = "logical",
        format = "character",
        db = "FAERSdbOrNULL"
    ),
    prototype = list(
        year = integer(),
        quarter = character(),
        data = NULL,
        meddra = NULL,
        deduplication = FALSE,
        standardization = FALSE,
        db = NULL
    )
)

################ validate method ########################
validate_faers <- function(object) {
    if (length(object@year) != length(object@quarter)) {
        return("the length of `@year` and `@quarter` must be the same")
    }
    if (!all(object@quarter %chin% FAERS_FILE_QUARTERS)) {
        return(sprintf(
            "`@quarter` must be values in %s",
            oxford_comma(FAERS_FILE_QUARTERS, final = "or")
        ))
    }
    if (anyDuplicated(faers_period(object))) {
        return("the period combined from `@year` and `@quarter` must be unique, you cannot import duplicated FAERS Quarterly Data file")
    }

    if (length(object@standardization) != 1L || is.na(object@standardization)) {
        return("@standardization must be a bool, `TRUE` or `FALSE`")
    }
    if (length(object@deduplication) != 1L || is.na(object@deduplication)) {
        return("@deduplication must be a bool, `TRUE` or `FALSE`")
    }
    if (!rlang::is_string(object@format, FAERS_FILE_FORMAT)) {
        return(sprintf(
            "`@format` must be a string of %s",
            oxford_comma(FAERS_FILE_FORMAT, final = "or")
        ))
    }
    TRUE
}

methods::setValidity("FAERS", validate_faers)

#' @export
#' @aliases FAERSascii
#' @rdname FAERS-class
methods::setClass(
    "FAERSascii",
    slots = list(data = "list", deletedCases = "character"),
    prototype = list(
        data = list(),
        deletedCases = character(),
        format = "ascii"
    ),
    contains = "FAERS"
)

methods::setValidity("FAERSascii", function(object) {
    data <- object@data
    if (!all(FAERS_ASCII_FILE_FIELDS %chin% names(data))) {
        return(sprintf(
            "`@data` must contain the all ascii fields, including %s", oxford_comma(FAERS_ASCII_FILE_FIELDS)
        ))
    }
    if (is.null(object@db)) {
        # memory mode: every field is a real data.table
        if (!all(vapply(data, data.table::is.data.table, logical(1L)))) {
            return("`@data` must be a list of `data.table`")
        }
    } else {
        # duckdb mode: every field is a FAERSdbTbl proxy backed by @db
        if (!all(vapply(data, methods::is, logical(1L), class2 = "FAERSdbTbl"))) {
            return("`@data` must be a list of `FAERSdbTbl` when `@db` is set")
        }
    }
    validate_faers(object)
})

#' @export
#' @aliases FAERSxml
#' @rdname FAERS-class
methods::setClass(
    "FAERSxml",
    slots = list(data = "DTOrNull", header = "list"),
    prototype = list(
        data = NULL,
        header = list(), format = "xml"
    ),
    contains = "FAERS"
)

################ show method ########################
#' @param object A [FAERS] object.
#' @importFrom methods show
#' @export
#' @method show FAERS
#' @rdname FAERS-class
methods::setMethod("show", "FAERS", function(object) {
    l <- nrow(faers_period(object))
    msg <- sprintf(
        "FAERS data from %s Quarterly %s file%s",
        l, object@format, if (l > 1L) "s" else ""
    )
    if (object@standardization && object@deduplication) {
        prefix <- "Standardized and De-duplicated"
    } else if (object@standardization) {
        prefix <- "Standardized"
    } else if (object@deduplication) {
        prefix <- "De-duplicated"
    } else {
        prefix <- NULL
    }
    if (!is.null(prefix)) msg <- paste(prefix, msg, sep = " ")
    cat(msg, sep = "\n")
    invisible(object)
})

#' @export
#' @method show FAERSascii
#' @rdname FAERS-class
methods::setMethod("show", "FAERSascii", function(object) {
    methods::callNextMethod(object)
    if (is.null(object@db)) {
        n_reports <- nrow(object@data$demo)
    } else {
        n_reports <- db_nrow(object@db@con, "demo")
    }
    if (object@deduplication) {
        msg <- sprintf(
            "  Total unique report%s: %s",
            if (n_reports > 1L) "s" else "", n_reports
        )
    } else {
        msg <- sprintf(
            "  Total report%s: %s (with duplicates)",
            if (n_reports > 1L) "s" else "", n_reports
        )
    }
    cat(msg, sep = "\n")
    invisible(object)
})

#######################################################
#' @param object A `FAERS` object.
#' @param ... Other arguments passed to specific methods.
#' @export
#' @aliases faers_data
#' @rdname FAERS-class
methods::setGeneric("faers_data", function(object, ...) {
    standardGeneric("faers_data")
})

#' @export
#' @method faers_data FAERS
#' @rdname FAERS-class
methods::setMethod("faers_data", "FAERS", function(object) {
    if (!is.null(object@db) && methods::is(object, "FAERSascii")) {
        return(materialize(object)@data)
    }
    object@data
})

#' @export
#' @aliases faers_year
#' @rdname FAERS-class
methods::setGeneric("faers_year", function(object) {
    standardGeneric("faers_year")
})

#' @export
#' @method faers_year FAERS
#' @rdname FAERS-class
methods::setMethod("faers_year", "FAERS", function(object) {
    object@year
})

##################################################################
#' @export
#' @aliases faers_quarter
#' @rdname FAERS-class
methods::setGeneric("faers_quarter", function(object) {
    standardGeneric("faers_quarter")
})

#' @export
#' @method faers_quarter FAERS
#' @aliases faers_quarter
#' @rdname FAERS-class
methods::setMethod("faers_quarter", "FAERS", function(object) {
    object@quarter
})

#################################################################
#' @export
#' @aliases faers_period
#' @rdname FAERS-class
methods::setGeneric("faers_period", function(object) {
    standardGeneric("faers_period")
})

#' @export
#' @method faers_period FAERS
#' @aliases faers_period
#' @rdname FAERS-class
methods::setMethod("faers_period", "FAERS", function(object) {
    data.table(year = object@year, quarter = object@quarter)
})

#' @export
#' @aliases faers_meddra
#' @rdname FAERS-class
methods::setGeneric("faers_meddra", function(object, ...) {
    standardGeneric("faers_meddra")
})

#' @param use A string, what meddra data to use, "hierarchy" or "smq". If
#' `NULL`, a [MedDRA] will be returned. Only used when `object` has been
#' standardized
#' @export
#' @method faers_meddra FAERS
#' @rdname FAERS-class
methods::setMethod("faers_meddra", "FAERS", function(object, use = NULL) {
    if (object@standardization) {
        out <- object@meddra
        if (!is.null(use)) {
            out <- methods::slot(out, match.arg(use, c("hierarchy", "smq")))
        }
    } else {
        out <- NULL
    }
    out
})

#################################################################
#' @export
#' @aliases faers_deleted_cases
#' @rdname FAERS-class
methods::setGeneric("faers_deleted_cases", function(object, ...) {
    standardGeneric("faers_deleted_cases")
})

#' @export
#' @method faers_deleted_cases FAERS
#' @aliases faers_field
#' @rdname FAERS-class
methods::setMethod("faers_deleted_cases", "FAERSascii", function(object) {
    object@deletedCases
})

#######################################################
#' @export
#' @aliases faers_header
#' @rdname FAERS-class
methods::setGeneric("faers_header", function(object) {
    standardGeneric("faers_header")
})

#' @export
#' @method faers_header FAERSxml
#' @aliases faers_header
#' @rdname FAERS-class
methods::setMethod("faers_header", "FAERSxml", function(object) {
    object@header
})
