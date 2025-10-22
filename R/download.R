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
    format <- match.arg(format, faers_file_format)
    assert_string(dir, empty_ok = FALSE)
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
faers_file_format <- c("ascii", "xml")
download_inform <- function(urls, file_paths, ...) {
    out <- file_paths
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
        
        old_timeout <- getOption("timeout")
        options(timeout = 300) 
        
        status <- utils::download.file(urls,
                                       destfile = file_paths, ...,
                                       method = "libcurl"
        )
        
        options(timeout = old_timeout)
        
        is_success <- status == 0L
        is_need_deleted <- !is_success & file.exists(file_paths)
        if (any(is_need_deleted)) file.remove(file_paths[is_need_deleted])
        if (!all(is_success)) {
            n_failed_files <- sum(!is_success)
            cli::cli_abort(c(
                "Cannot download {.val {n_failed_files}} file{?s}",
                "i" = "url{?s}: {.url {urls[!is_success]}}",
                "!" = paste(
                    "status {cli::qty(n_failed_files)} code{?s}:",
                    "{.val {status[!is_success]}}"
                )
            ))
        }
    }
    out
}

#' @param status A data frame returned by [multi_download][curl::multi_download]
#' @noRd
is_download_success <- function(status, successful_code = c(200L, 206L, 416L)) {
    !is.na(status$success) &
        status$success &
        (status$status_code %in% successful_code)
}

base_download_inform <-  function(urls, file_paths, ...) {
    out <- file_paths
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
        
        if (!requireNamespace("httr", quietly = TRUE)) {
            stop("please install httr package: install.packages('httr')")
        }
        
        user_agent <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
        
        download_with_retry <- function(url, destfile, max_attempts = 3) {
            for (attempt in 1:max_attempts) {
                tryCatch({
                    if (attempt > 1) Sys.sleep(2)
                    
                    cli::cli_inform("try downloading {attempt}/{max_attempts}")
                    
                    response <- httr::GET(
                        url,
                        httr::user_agent(user_agent),
                        httr::add_headers(
                            Accept = "application/octet-stream",
                            Referer = "https://www.fda.gov/"
                        ),
                        httr::timeout(30)
                    )
                    
                    if (httr::status_code(response) != 200) {
                        stop("HTTP Wrong: ", httr::status_code(response))
                    }
                    
                    writeBin(httr::content(response, "raw"), destfile)
                    
                    if (file.exists(destfile) && file.size(destfile) > 0) {
                        cli::cli_inform("✓ success: {.file {basename(destfile)}}")
                        return(0L)  
                    } else {
                        stop("The downloaded file is empty or does not exist.")
                    }
                    
                }, error = function(e) {
                    cli::cli_warn("try downloading {attempt} fail: {e$message}")
                    
                    if (file.exists(destfile)) file.remove(destfile)
                    if (attempt == max_attempts) {
                        stop("All download attempts have failed.: ", e$message)
                    }
                })
            }
            return(1L)  
        }
        
        status <- mapply(download_with_retry, urls, file_paths)
        
        is_success <- status == 0L
        is_need_deleted <- !is_success & file.exists(file_paths)
        if (any(is_need_deleted)) file.remove(file_paths[is_need_deleted])
        
        if (!all(is_success)) {
            n_failed_files <- sum(!is_success)
            cli::cli_abort(c(
                "can't download {.val {n_failed_files}} file",
                "i" = "URL: {.url {urls[!is_success]}}",
                "x" = "Please check.:",
                " " = "- Network connection",
                " " = "- Whether the URL is valid",
                " " = "- Website access permissions"
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
