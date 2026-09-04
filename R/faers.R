#' Download and parse FAERS Quarterly Data files
#' @inheritParams faers_download
#' @inheritParams faers_parse
#' @param handle_opts Extra handle options passed to each request
#' [new_handle][curl::new_handle].
#' @return A [FAERSxml] or [FAERSascii] object.
#' @examples
#' # you must change `dir`, as the file included in the package is sampled
#' data <- faers(2004, "q1",
#'     dir = system.file("extdata", package = "faers"),
#'     compress_dir = tempdir()
#' )
#' @param database A string, one of `"memory"` (default) or `"duckdb"`.  With
#' `"memory"`, all quarters are held in memory (original behavior).  With
#' `"duckdb"`, quarterly data is written into a DuckDB database so it does not
#' stay resident in memory; you keep the same analysis pipeline
#' (standardize → dedup → counts → phv_signal) unchanged.
#' @param db_path A string, the DuckDB file path used when `database =
#' "duckdb"`.  Defaults to `":memory:"` (a transient in-memory database).
#' @export
faers <- function(years, quarters, format = NULL, dir = getwd(), compress_dir = dir, handle_opts = list(),
                  database = c("memory", "duckdb"), db_path = NULL) {
    database <- match.arg(database)
    format <- match.arg(format, FAERS_FILE_FORMAT)
    yq <- recycle_scalar(years = years, quarters = quarters)
    data.table::setDT(yq)
    yq <- unique(yq)
    faers_files <- do.call(faers_download, c(
        yq, list(format = format, dir = dir), handle_opts
    ))
    bar_id <- cli::cli_progress_bar(
        "Parsing FAER",
        type = "iterator", total = length(faers_files),
        format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}",
        format_done = sprintf(
            "Parsing {.val {cli::pb_total}} %s Quarterly Data file{?s} in {cli::pb_elapsed}",
            format
        ),
        clear = FALSE
    )
    # For `database = "duckdb"`, when more than one quarter is parsed into a
    # shared on-disk file, each quarter needs its own database file; otherwise
    # the per-quarter `db_ingest(overwrite=TRUE)` calls would clobber each
    # other.  (:memory: needs no handling — each quarter gets its own DB, and
    # combine_faers_db merges them into a fresh file DB.)
    if (database == "duckdb" && length(yq$years) > 1L &&
        !is.null(db_path) && !identical(db_path, ":memory:")) {
        per_quarter_db <- file.path(
            tempdir(),
            sprintf("faers_quarter_%03d_%s.duckdb", seq_len(nrow(yq)), digest_short(db_path))
        )
    } else {
        per_quarter_db <- rep(list(db_path), nrow(yq))
    }
    out <- .mapply(
        function(path, year, quarter, format, compress_dir, index) {
            out <- faers_parse(
                path = path,
                format = format, year = year, quarter = quarter,
                compress_dir = compress_dir,
                database = database,
                db_path = per_quarter_db[[index]]
            )
            cli::cli_progress_update(id = bar_id)
            out
        },
        list(
            path = faers_files, year = yq$years, quarter = yq$quarters, index = seq_len(nrow(yq))
        ),
        MoreArgs = list(format = format, compress_dir = compress_dir)
    )
    l <- length(out)
    if (l == 1L) {
        return(out[[1L]])
    }
    cli::cli_alert("Combining all {l} {.cls FAERS} Data{?s}")
    combine_faers(out)
}

recycle_scalar <- function(..., length = NULL, args = NULL) {
    args <- args %||% unlist(lapply(substitute(...()), as.character),
        use.names = FALSE
    )
    lst <- list(...)
    l <- lengths(lst)
    expected_len <- length %||% max(l)
    if (!all(l == 1L | l == expected_len)) {
        cli::cli_abort(c(
            "{.arg {args}} must have compatible sizes",
            i = "Only values of size one are recycled."
        ))
    }
    lapply(lst, rep_len, length.out = expected_len)
}
