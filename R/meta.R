#' List of FAERS data
#'
#' The function lists the metadata for the FAERS databases currently
#' available to download.
#'
#' @param force A boolean value. If set to `TRUE`, it indicates the retrieval of
#' information about all records' metadata in the FAERS Quarterly Data Extract
#' Files Site, bypassing the cache.
#' @param internal A boolean value. It determines whether to use the internal
#' data associated with the package.
#' @return A [data.table][data.table::data.table] reporting years, period,
#' quarter, and file urls and file sizes.
#' @examples
#' # Queries the FDA website (falls back to bundled metadata without internet):
#' faers_meta(internal = !curl::has_internet())
#' @seealso
#' <https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html>
#' @export
faers_meta <- function(force = FALSE, internal = !curl::has_internet()) {
    assert_bool(force)
    assert_bool(internal)
    if (force || is.null(out <- faers_meta_cache_read(internal = internal))) {
        out <- faers_meta_parse()
        faers_meta_cache_save(out)
    }
    out
}

faers_meta_cache_read <- function(internal = FALSE) {
    if (internal) {
        out <- load_data("faers_meta_data")
        msg <- "Using internal FAERS metadata"
    } else if (file.exists(file <- faers_meta_cache_file())) {
        out <- readRDS(file)
        msg <- "Using FAERS metadata from cached {.file {file}}"
        # save the data in the cached environment for the usage of next time
        # like faers_available()
    } else {
        return(NULL)
    }
    cli::cli_inform(c(">" = msg, " " = "Snapshot time: {out$date}"))
    out$data
}

faers_meta_cache_save <- function(data) {
    file <- faers_meta_cache_file()
    cli::cli_inform(
        c(">" = "Writing FAERS metadata into cache file {.file {file}}")
    )
    saveRDS(list(data = data, date = Sys.time()), file = file)
}

faers_meta_cache_file <- function(dir = faers_cache_dir("metadata")) {
    file.path(dir, "faers_meta_data.rds")
}

faers_meta_parse <- function(call = rlang::caller_env()) {
    assert_internet(call = call)
    url <- sprintf(
        "%s/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html",
        fda_host("fis")
    )
    cli::cli_inform(c(">" = "Reading html: {.url {url}}"))
    html <- xml2::read_html(url)
    table_xml_list <- rvest::html_elements(html, ".panel.panel-default")
    table_list <- lapply(table_xml_list, parse_xml_table)
    out <- data.table::rbindlist(table_list)
    out[, quarter := period2quarter(period)] # nolint
    data.table::setcolorder(out, c(
        "year", "quarter", "period",
        "urls_ascii", "file_size_ascii", "urls_xml", "file_size_xml"
    ))
    data.table::setnames(out, c(
        "year", "quarter", "period",
        "ascii_urls", "ascii_file_size", "xml_urls", "xml_file_size"
    ))
    data.table::setorderv(out, c("year", "quarter"), order = c(-1L, -1L))[]
}

utils::globalVariables(c("period", "quarter"))

parse_xml_table <- function(year_xml) {
    year <- rvest::html_text2(rvest::html_element(year_xml, ".panel-title"))
    if (!str_detect(year, "^20\\d+$")) {
        return(NULL)
    }
    file_table <- lapply(rvest::html_elements(year_xml, "tbody > tr"), function(quarter_xml) {
        period <- rvest::html_text2(rvest::html_elements(quarter_xml, "p"))
        file_xmls <- rvest::html_elements(quarter_xml, "td a")
        files <- rvest::html_text2(file_xmls)
        out <- data.table(
            period = str_remove(period, "\\s*\\d*(\\s*posted.*\\s*)?$"),
            urls = rvest::html_attr(file_xmls, "href"),
            file_format = tolower(str_extract(files, "XML|ASCII")),
            file_size = str_extract(files, "\\d[\\d.]*(MB)")
        )
        data.table::dcast(out, period ~ file_format,
            value.var = c("urls", "file_size")
        )
    })
    out_table <- data.table::rbindlist(file_table)
    out_table[, year := as.integer(year)] # nolint
}

period2quarter <- function(x) {
    data.table::fcase(
        str_detect(x, "^October"), "q4",
        str_detect(x, "^July"), "q3",
        str_detect(x, "^April"), "q2",
        str_detect(x, "^January"), "q1"
    )
}
