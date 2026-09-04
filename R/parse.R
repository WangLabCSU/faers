#' Parse FAERS Quarterly Data
#' @param path A string specifies the path of FAERS Quarterly Data. You can pass
#' the FAERS zip file directly, In this way, all files in the zip file will be
#' extracted in `compress_dir`. Or, you can also uncompressed youself, and
#' passed the directory contained the uncompressed files.
#' @inheritParams faers_download
#' @param year Year of the FAERS Quarterly Data. Coerced into integer, if
#' `NULL`, this will be extracted from path.
#' @param quarter String specifies quarter of the FAERS data, if `NULL`, this
#' will be extracted from path.
#' @param compress_dir A string specifies the directory to extract files to. It
#' will be created if necessary.
#' @return A [FAERSxml] or [FAERSascii] object.
#' @examples
#' # the files included in the package are sampled
#' data <- faers_parse(
#'     system.file("extdata", "aers_ascii_2004q1.zip", package = "faers"),
#'     compress_dir = tempdir()
#' )
#' @param database A string, one of `"memory"` (default) or `"duckdb"`.  With
#' `"memory"`, data is held in memory exactly as before.  With `"duckdb"`, the
#' quarter data is written into a DuckDB file so it does not stay resident in
#' memory; see [faers] for details.
#' @param db_path A string, the DuckDB file path used when `database =
#' "duckdb"`.  Defaults to `":memory:"` (a transient in-memory database).
#' @export
faers_parse <- function(path, format = NULL, year = NULL, quarter = NULL, compress_dir = getwd(),
                        database = c("memory", "duckdb"), db_path = NULL) {
    database <- match.arg(database)
    assert_string(path, allow_empty = FALSE)
    if (is.null(format)) {
        format <- str_extract(basename(path), "xml|ascii", ignore.case = TRUE)
        if (!any(format == FAERS_FILE_FORMAT)) {
            cli::cli_abort(c(
                "Cannot parse file format from {.arg path}",
                i = "Try to set {.arg format} manually"
            ))
        }
    } else {
        format <- match.arg(format, FAERS_FILE_FORMAT)
    }
    year <- year %||% str_extract(path, "20\\d+(?=q[1-4])")
    year <- as.integer(year)
    quarter <- quarter %||% str_extract(path, "(?<=20\\d{2})(q[1-4])")
    quarter <- as.character(quarter)
    path <- dir_or_unzip(path,
        compress_dir = compress_dir,
        pattern = "20\\d{2}q[1-4]\\.zip$",
        none_msg = c(
            "Only compressed zip files from FAERS Quarterly Data can work",
            i = "with pattern: \"20\\d{{2}}q[1-4]\\.zip\""
        )
    )
    if (database == "duckdb") {
        return(parse_to_db(path, year, quarter, format = format, db_path = db_path))
    }
    switch(format,
        xml = parse_xml(path, year, quarter),
        ascii = parse_ascii(path, year, quarter)
    )
}

# parse ascii files -----------------------------------
parse_ascii <- function(path, year, quarter) {
    files <- locate_files(locate_dir(path, "^ascii$"), "\\.txt$")
    # for 2018q1 demo file, there exists a suffix "_new"
    fields <- tolower(str_remove(
        basename(files), "\\d+q\\d(_new)?\\.txt$",
        ignore.case = TRUE
    ))
    idx <- match(FAERS_ASCII_FILE_FIELDS, fields)
    if (anyNA(idx)) {
        cli::cli_abort(sprintf(
            "Cannot find %s",
            oxford_comma(style_file(FAERS_ASCII_FILE_FIELDS[is.na(idx)]))
        ))
    }
    files <- files[idx]
    fields <- fields[idx]
    deleted_cases <- read_ascii_deleted_cases(path, year, quarter)
    data_list <- .mapply(function(file, field) {
        if (any(year == ascii_file_need_safe_parsing$year &
            quarter == ascii_file_need_safe_parsing$quarter &
            field == ascii_file_need_safe_parsing$field)) {
            out <- safely_read_ascii(file, year, quarter, field = field)
        } else {
            out <- read_ascii(file, verbose = FALSE)
        }
        unify_ascii(out, field = field, year = year, quarter = quarter)
        out[!is.na(primaryid)]
    }, list(file = files, field = fields), NULL)
    data.table::setattr(data_list, "names", fields)
    methods::new("FAERSascii",
        data = data_list,
        year = year, quarter = quarter,
        deletedCases = deleted_cases
    )
}

read_ascii_deleted_cases <- function(path, year, quarter) {
    # As of 2019 Quarter one there is a new text file that lists deleted files
    if (!is_before_period(year, quarter, 2018L, "q4")) {
        deleted_cases_files <- locate_files(locate_dir(path, "^deleted$"))
        deleted_cases <- lapply(deleted_cases_files, function(file) {
            data.table::fread(
                file = file,
                header = FALSE,
                na.strings = na_string,
                blank.lines.skip = TRUE,
                keepLeadingZeros = TRUE
            )[[1L]]
        })
        unique(as.character(unlist(deleted_cases, use.names = FALSE)))
    } else {
        character()
    }
}

read_ascii <- function(file, ...) {
    out <- data.table::fread(
        file = file,
        sep = "$", quote = "", fill = TRUE,
        blank.lines.skip = TRUE,
        na.strings = na_string,
        keepLeadingZeros = TRUE,
        ...
    )
    # the last columns often messed by the presence OF '$'.
    # AS PART OF THE ROWTERMINATOR IN DATA ROWS BUT NOT IN THE HEADER ROW
    # (LAERS)
    last_col <- names(out)[ncol(out)]
    if (str_detect(last_col, "^V\\d+$")) {
        out[, (last_col) := lapply(.SD, function(x) {
            if (all(is.na(x))) NULL else x
        }), .SDcols = last_col]
    }
    out
}


# data.table will stop early for 4 files
# leave a lot of rows not read in
# this is mainly due to the presence of the collapsed lines (two, or more
# lines collapsed as one line)
# aers_ascii_2012q1/ascii/DEMO12Q1.TXT
# aers_ascii_2011q4/ascii/DRUG11Q4.TXT
# aers_ascii_2011q3/ascii/DRUG11Q3.TXT
# aers_ascii_2011q2/ascii/DRUG11Q2.TXT
ascii_file_need_safe_parsing <- data.table(
    year = c(2011L, 2011L, 2011L, 2012L),
    quarter = c("q2", "q3", "q4", "q1"),
    field = c("drug", "drug", "drug", "demo")
)

safely_read_ascii <- function(file, year, quarter, field) {
    file_text <- brio::read_lines(file)
    if (year == 2011L && quarter == "q2" && field == "drug") {
        file_text[322967L] <- str_replace(
            # line: 322967
            file_text[322967L],
            "7475791$1016572490$SS$DOXORUBICIN",
            "\n7475791$1016572490$SS$DOXORUBICIN",
            fixed = TRUE
        )
    }
    if (year == 2011L && quarter == "q3" && field == "drug") {
        file_text[247896L] <- str_replace(
            # line: 247896
            file_text[247896L],
            "7652730$1017255397$SS$BEVACIZUMAB",
            "\n7652730$1017255397$SS$BEVACIZUMAB",
            fixed = TRUE
        )
    }
    if (year == 2011L && quarter == "q4" && field == "drug") {
        file_text[446738L] <- str_replace(
            file_text[446738L],
            "7941354$1018188213$SS$MEMANTINE HYDROCHLORIDE",
            "\n7941354$1018188213$SS$MEMANTINE HYDROCHLORIDE",
            fixed = TRUE
        )
    }
    # for 2012q1 demo data:
    # JP-CUBIST-$E2B0000000182, the "$" should be removed
    if (year == 2012L && quarter == "q1" && field == "demo") {
        file_text[105917L] <- str_replace(
            file_text[105917L],
            "JP-CUBIST-$E2B0000000182",
            "JP-CUBIST-E2B0000000182",
            fixed = TRUE
        )
    }

    # n_seps <- str_count(file_text, "$", fixed = TRUE)
    # collapsed_lines <- floor(n_seps / n_seps[2L]) > 1L
    # if (any(collapsed_lines)) {
    #     cli::cli_warn(c(
    #         "omiting {sum(collapsed_lines)} collapsed line{?s}",
    #         i = "line number: {which(collapsed_lines)}"
    #     ))
    #     file_text <- file_text[!collapsed_lines]
    # }
    # file_text <- str_remove(file_text[!collapsed_lines], "\\$$")
    read_text(file_text,
        sep = "$", quote = "", fill = TRUE,
        blank.lines.skip = TRUE,
        integer64 = "double"
    )
}

read_text <- function(text, ...) {
    if (!length(text)) {
        return(data.table())
    }
    file <- tempfile()
    data.table::fwrite(list(text),
        file = file,
        quote = FALSE,
        na = "NA",
        col.names = FALSE,
        logical01 = FALSE,
        showProgress = FALSE,
        compress = "none",
        verbose = FALSE
    )
    # brio::write_lines(text, file)
    on.exit(file.remove(file))
    data.table::fread(
        file = file, ...,
        na.strings = na_string,
        showProgress = FALSE,
        keepLeadingZeros = TRUE
    )
}
na_string <- c("NA", "")

# fix for 2011q4 drug file, there exists a \177 field in line 729342
# which induce fread stop early, so I just don't use this function.
read_lines <- function(file) {
    data.table::fread(
        file = file, sep = "", header = FALSE,
        blank.lines.skip = TRUE,
        colClasses = "character",
        showProgress = FALSE,
        keepLeadingZeros = TRUE
    )[[1L]]
}

# parse xml ----------------------------------------------
parse_xml <- function(path, year, quarter) {
    file <- locate_files(locate_dir(path, "^xml$"), "\\.xml$")
    xml_doc <- xml2::read_xml(file)
    full_content <- xml2::xml_contents(xml_doc)
    header <- xml2::as_list(full_content[[1L]])
    reports_nodesets <- full_content[-1L]
    l <- length(reports_nodesets)
    bar_id <- cli::cli_progress_bar("Parsing reports",
        type = "iterator", total = l,
        format = "{cli::pb_spin} Parsing reports | {cli::pb_current}/{cli::pb_total}",
        format_done = "Parsing {.val {cli::pb_total}} report{?s} in {cli::pb_elapsed}",
        clear = FALSE
    )
    reports_list <- lapply(seq_len(l), function(i) {
        report <- xml2::as_list(reports_nodesets[[i]])
        report <- make_length_one_list(list_flatten(report))
        for (item in c("patient.drug", "patient.reaction")) {
            idx <- names(report) == item
            item_list <- report[idx]
            report <- report[!idx]
            lapply(item_list, function(x) data.table::setDT(x[[1L]]))
            item_dt <- data.table::rbindlist(
                unlist(item_list, recursive = FALSE, use.names = FALSE),
                use.names = TRUE, fill = TRUE
            )
            report[[item]] <- list(simplify_list_cols(item_dt))
        }
        cli::cli_progress_update(id = bar_id)
        data.table::setDT(report)
    })
    data <- data.table::rbindlist(reports_list,
        use.names = TRUE, fill = TRUE
    )
    simplify_list_cols(data)
    data.table::setnames(data, tolower)
    methods::new("FAERSxml",
        data = data, header = header,
        year = year, quarter = quarter
    )
}

list_flatten <- function(x, use.names = TRUE, sep = ".") {
    lst <- unlist(unname(x), recursive = FALSE)
    if (!use.names) {
        return(lst)
    }
    top_names <- rep(rlang::names2(x), times = lengths(x))
    is_top_names_good <- top_names != ""
    cur_names <- rlang::names2(lst)
    is_cur_names_good <- cur_names != ""
    new_names <- data.table::fcase(
        is_top_names_good & is_cur_names_good,
        paste(top_names, cur_names, sep = sep),
        is_top_names_good, top_names,
        is_cur_names_good, cur_names,
        default = NA_character_
    )
    data.table::setattr(lst, "names", new_names)
}

make_length_one_list <- function(lst) {
    lapply(lst, function(x) {
        if (length(x) == 1L) x else list(x)
    })
}

simplify_list_cols <- function(dt) {
    dt_names <- names(dt)
    dt[, (dt_names) := lapply(.SD, function(x) { # nolint
        x[lengths(x) == 0L] <- NA
        x[lengths(x) == 1L] <- unlist(x[lengths(x) == 1L])
        x
    }), .SDcols = dt_names]
    flatten_cols <- dt_names[
        vapply(dt, function(x) all(lengths(x) == 1L), logical(1L))
    ]
    dt[, (flatten_cols) := lapply(.SD, unlist, recursive = TRUE), # nolint
        .SDcols = flatten_cols
    ]
}
