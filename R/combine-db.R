# Database-mode combination of multiple quarterly FAERS objects.
#
# combine_faers_db takes several db-backed FAERSascii objects (each holding a
# distinct quarter) and produces the @data proxies + @db store for a single
# combined object whose DuckDB tables hold every quarter's rows.  This mirrors
# combine_faers_ascii_data (rbindlist over fields) but works field-by-field:
# each field is read from every source, bound together, and written to a fresh
# database before the next field is read.  Only one field is resident at a
# time, so peak memory stays bounded even when combining many quarters.

combine_faers_db <- function(x) {
    # All elements are guaranteed (by combine_faers dispatch) to be db-backed.
    .path <- db_new_path()
    con <- db_connect(.path)

    .first <- x[[1L]]
    fields <- names(.first@db@tables)
    data <- list()
    for (field in fields) {
        parts <- lapply(x, function(obj) db_collect(obj@db@con, field))
        combined <- data.table::rbindlist(parts, fill = TRUE, use.names = TRUE)
        rm(parts); gc()
        db_ingest(con, field, combined)
        rm(combined); gc()
        data[[field]] <- methods::new("FAERSdbTbl", con = con, table = db_field_table(field))
    }

    db_store <- methods::new(
        "FAERSdb",
        con = con,
        path = .path,
        version = tryCatch(
            as.character(utils::packageVersion("duckdb")),
            error = function(e) NA_character_
        ),
        tables = db_field_tables()
    )
    list(data = data, db = db_store)
}
