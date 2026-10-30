#' Check if FAERS year
#'
#' This function check if data for the years and quarters selected are
#' available at FAERS to be downloaded.
#'
#' @param years An atomic integer indicates years for which data are required.
#' @param quarters An atomic character, only "q1", "q2", "q3", and "q4" are
#' allowed.
#' @inheritParams faers_meta
#' @return A logical indicates FAERS can have data for the `years` and
#' `quarters` required?
#' @examples
#' # Queries the FDA website (falls back to bundled metadata without internet):
#' faers_available(c(2011, 2023), c("q1", "q2"), internal = !curl::has_internet())
#' @export
faers_available <- function(years, quarters, force = FALSE, internal = FALSE) {
    assert_inclusive(quarters, FAERS_FILE_QUARTERS)
    metadata <- faers_meta(force = force, internal = internal)
    years %in% metadata$year & quarters %chin% metadata$quarter
}
