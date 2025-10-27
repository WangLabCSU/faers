#' Download FAERS data
#'
#' This function downloads the FAERS data for selected years and quarters.
#'
#' @inheritParams faers_available
#' @param format File format to used, only "ascii" and "xml" are availabe.
#'  Default: "ascii".
#' @param dir The destination directory for any downloads. Defaults to
#'  current working dir.
#' @param ... Extra handle options passed to each request
#'  [multi_download][curl::multi_download].
#' @return An atomic character for the path of downloaded files.
#' @examples
#' # you must change `dir`, as the file included in the package is sampled
#' # in this way, the file will downloaded from FAERS
#' faers_download(
#'     year = 2004, quarter = "q1",
#'     dir = system.file("extdata", package = "faers")
#' )
#' @export
faers_download <- function(years, quarters, format = NULL, dir = getwd(), ...) {
    format <- match.arg(format, FAERS_FILE_FORMAT)
    assert_string(dir, allow_empty = FALSE)
    if (format == "xml") {
        # only faers database has xml data files
        is_aers_pairs <- is_from_laers(years, quarters)
        if (any(is_aers_pairs)) {
            aers_pairs <- paste0(years, quarters)[is_aers_pairs] # nolint
            cli::cli_abort(c(
                "Only FAERS (from 2012q4) has {.field xml} files",
                x = "Legacy AERS (before 2012q3) pair{?s}: {.val {aers_pairs}}"
            ))
        }
    }
    urls <- build_faers_url(format, years, quarters)
    dest_files <- file.path(dir_create2(dir), basename(urls))
    download_inform(urls, dest_files, ...)
}

#' Download utils function with good message.
#' @return A character path if downloading successed, otherwise, stop with error
#'   message.
#' @noRd
download_inform <- function(urls, file_paths, ...) {
    ans <- file_paths
    if (any(is_existed <- file.exists(file_paths))) {
        cli::cli_inform(paste(
            "Finding {.val {sum(is_existed)}} file{?s} already",
            "downloaded: {.file {basename(file_paths[is_existed])}}"
        ))
        urls <- urls[!is_existed]
        file_paths <- file_paths[!is_existed]
    }
    if (l <- length(urls)) {
        assert_internet()
        if (l == 1L) {
            cli::cli_inform("Downloading 1 file from: {.url {urls}}")
        } else {
            cli::cli_inform("Downloading {.val {l}} files")
        }
        
        headers <- get_complete_headers()
        
        arg_list <- c(
            list(
                urls = urls, destfiles = file_paths, resume = FALSE,
                progress = interactive(), multi_timeout = Inf
            ),
            list(
                useragent = headers$`User-Agent`,
                httpheader = c(
                    Accept = headers$Accept,
                    `Accept-Language` = headers$`Accept-Language`,
                    `Accept-Encoding` = headers$`Accept-Encoding`,
                    Connection = headers$Connection,
                    Referer = headers$Referer
                )
            ),
            rlang::list2(...)
        )
        status <- do.call(curl::multi_download, arg_list)
        is_success <- is_download_success(status)
        is_need_deleted <- !is_success & file.exists(file_paths)
        if (any(is_need_deleted)) file.remove(file_paths[is_need_deleted])
        if (!all(is_success)) {
            n_failed_files <- sum(!is_success) # nolint
            cli::cli_abort(c(
                "Cannot download {.val {n_failed_files}} file{?s}",
                "i" = "url{?s}: {.url {urls[!is_success]}}",
                "!" = paste(
                    "status {cli::qty(n_failed_files)} code{?s}:",
                    "{.val {status$status_code[!is_success]}}"
                ),
                x = paste(
                    "error {cli::qty(n_failed_files)} message{?s}:",
                    "{.val {status$error[!is_success]}}"
                )
            ))
        }
    }
    ans
}

#' @param status A data frame returned by [multi_download][curl::multi_download]
#' @noRd
is_download_success <- function(status, successful_code = c(200L, 206L, 416L)) {
    !is.na(status$success) &
        status$success &
        (status$status_code %in% successful_code)
}

base_download_inform <- function(urls, file_paths, ...) {
    out <- file_paths
    if (any(is_existed <- file.exists(file_paths))) {
        cli::cli_inform(paste(
            "Finding {.val {sum(is_existed)}} file{?s} already",
            "downloaded: {.file {basename(file_paths[is_existed])}}"
        )) # nolint
        urls <- urls[!is_existed]
        file_paths <- file_paths[!is_existed]
    }
    if (l <- length(urls)) {
        assert_internet()
        if (l == 1L) {
            cli::cli_inform("Downloading file from: {.url {urls}}")
        } else {
            cli::cli_inform("Downloading {.val {l}} files")
        }
        

        headers <- get_complete_headers()
        
        arg_list <- list(
            urls = urls, destfiles = file_paths, resume = FALSE,
            progress = interactive(), multi_timeout = Inf,
            useragent = headers$`User-Agent`,
            httpheader = c(
                Accept = headers$Accept,
                `Accept-Language` = headers$`Accept-Language`,
                `Accept-Encoding` = headers$`Accept-Encoding`,
                Connection = headers$Connection,
                Referer = headers$Referer
            )
        )
        
        status <- do.call(curl::multi_download, arg_list)
        is_success <- is_download_success(status)
        is_need_deleted <- !is_success & file.exists(file_paths)
        if (any(is_need_deleted)) file.remove(file_paths[is_need_deleted])
        if (!all(is_success)) {
            n_failed_files <- sum(!is_success) # nolint
            cli::cli_abort(c(
                "Cannot download {.val {n_failed_files}} file{?s}",
                "i" = "url{?s}: {.url {urls[!is_success]}}",
                "!" = paste(
                    "status {cli::qty(n_failed_files)} code{?s}:",
                    "{.val {status$status_code[!is_success]}}"
                )
            ))
        }
    }
    out
}

build_faers_url <- function(type, years, quarters) {
    laers_period <- is_from_laers(years, quarters)
    sprintf(
        "%s/content/Exports/%s_%s_%s%s.zip",
        fda_host("fis"),
        ifelse(laers_period, "aers", "faers"),
        ifelse(type == "ascii" | !laers_period, type, "sgml"),
        years, quarters
    )
}

#' Complete Request Header Module
#' @return 
#' @noRd
get_complete_headers <- function() {
    user_agents <- c(
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36",
        "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/14.1.1 Safari/605.1.15",
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:89.0) Gecko/20100101 Firefox/89.0",
        "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/92.0.4515.107 Safari/537.36",
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Edge/91.0.864.59"
    )
    

    referers <- c(
        "https://www.fda.gov/drugs/surveillance/fda-adverse-event-reporting-system-faers",
        "https://www.fda.gov/drugs/drug-approvals-and-databases/fda-adverse-event-reporting-system-faers",
        "https://www.fda.gov/drugs",
        "https://www.fda.gov",
        "https://google.com"
    )
    
    accept_encodings <- c(
        "gzip, deflate, br",
        "gzip, deflate",
        "gzip"
    )
    

    accept_languages <- c(
        "en-US,en;q=0.9",
        "en-US,en;q=0.8",
        "en-GB,en;q=0.9,en-US;q=0.8"
    )
    
    list(
        `User-Agent` = sample(user_agents, 1),
        Accept = "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
        `Accept-Language` = sample(accept_languages, 1),
        `Accept-Encoding` = sample(accept_encodings, 1),
        Connection = "keep-alive",
        Referer = sample(referers, 1)
    )
}
