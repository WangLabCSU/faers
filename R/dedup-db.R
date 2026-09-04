# Database-mode deduplication.
#
# Per design doc §3.6 (phase 1), dedup is a materializing transform: the data
# needed to *determine* which reports to keep (demo + the alignment columns of
# drug/indi/ther/reac) is pulled into memory, run through the *same*
# `dedup_faers_ascii`, and the resulting keep-set (year, quarter, primaryid)
# is used to trim every field's table in DuckDB.  Reusing `dedup_faers_ascii`
# verbatim guarantees the keep-set — and therefore the final data — is
# identical to the memory path.

dedup_db <- function(object, remove_deleted_cases = TRUE) {
    if (isTRUE(remove_deleted_cases)) {
        deleted_cases <- faers_deleted_cases(object)
        if (!length(deleted_cases)) deleted_cases <- NULL
    } else {
        deleted_cases <- NULL
    }

    # Pull the field data needed to run dedup into memory (as a list named as
    # `dedup_faers_ascii` expects). This is the same shape as faers_mget.
    data <- lapply(c("demo", "drug", "indi", "ther", "reac"), function(field) {
        db_collect(object@db@con, field)
    })
    names(data) <- c("demo", "drug", "indi", "ther", "reac")

    keep <- dedup_faers_ascii(data, deleted_cases = deleted_cases)
    rm(data); gc()

    # Trim every field table to the keep-set.
    db_keep_keys(object, keep)

    object@deduplication <- TRUE
    object
}

# Keep only rows whose (year, quarter, primaryid) is in `keep` (a data.table
# of those three columns). Rewrites each field's table in DuckDB.
db_keep_keys <- function(object, keep) {
    con <- object@db@con
    keep_df <- as.data.frame(keep)
    names(keep_df) <- c("year", "quarter", "primaryid")
    DBI::dbWriteTable(con, "__faers_keep_keys__", keep_df,
        overwrite = TRUE, temporary = TRUE
    )
    for (field in names(object@db@tables)) {
        table_name <- object@db@tables[[field]]
        DBI::dbExecute(
            con,
            sprintf(
                "DELETE FROM %s WHERE (year, quarter, primaryid) NOT IN
                     (SELECT year, quarter, primaryid FROM __faers_keep_keys__)",
                table_name
            )
        )
    }
    invisible(object)
}
