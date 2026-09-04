# Database-backed parsing: write FAERS quarter data into DuckDB instead of
# keeping all quarters resident in memory.
#
# v1 approach (per design doc §3.3, "read once, write to DB, release"):
#   * reuse the existing, well-tested parse_ascii / parse_xml to build the
#     unified data.tables for ONE quarter,
#   * immediately write them into the DuckDB file,
#   * release that quarter's in-memory copy.
# The headline memory win at this stage is architectural: instead of all
# quarters staying resident for the whole analysis (as `faers_combine` does
# with rbindlist), each quarter is parsed, written to disk, and freed.
#
# The resulting FAERSascii object carries an OPEN connection in @db; fields
# are materialized on demand (see methods-FAERS.R / counts / phv).

# Parse a single (unzipped) FAERS quarter and write it to DuckDB.
# Returns a FAERSascii with `@db` set.
parse_to_db <- function(path, year, quarter, format = c("ascii", "xml"),
                        db_path = NULL) {
    format <- match.arg(format, FAERS_FILE_FORMAT)
    db_require()
    con <- db_connect(db_path)
    path <- dir_or_unzip(path,
        compress_dir = getwd(),
        pattern = "20\\d{2}q[1-4]\\.zip$",
        none_msg = c(
            "Only compressed zip files from FAERS Quarterly Data can work",
            i = "with pattern: \"20\\d{{2}}q[1-4]\\.zip\""
        )
    )
    if (format == "ascii") {
        # Reuse the existing unified parse for ONE quarter.
        obj <- parse_ascii(path, year, quarter)
        fields <- names(obj@data)
        deleted_cases <- obj@deletedCases
        for (field in fields) {
            db_ingest(con, field, obj@data[[field]])
        }
    } else {
        obj <- parse_xml(path, year, quarter)
        # xml objects carry a single combined data table.
        fields <- "combined"
        deleted_cases <- character()
        db_ingest(con, "combined", obj@data)
    }

    tables <- stats::setNames(
        vapply(fields, db_field_table, character(1L)),
        fields
    )
    db_store <- methods::new(
        "FAERSdb",
        con = con,
        path = db_normalize_path(db_path),
        version = tryCatch(
            as.character(utils::packageVersion("duckdb")),
            error = function(e) NA_character_
        ),
        tables = tables
    )

    methods::new(
        "FAERSascii",
        data = stats::setNames(lapply(fields, function(f) {
            methods::new("FAERSdbTbl", con = con, table = db_field_table(f))
        }), fields),
        year = year, quarter = quarter,
        deletedCases = deleted_cases,
        standardization = FALSE, deduplication = FALSE,
        format = format,
        db = db_store
    )
}
