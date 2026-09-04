# Database-mode faers_counts: push the distinct-primaryid count down to SQL.
#
# The in-memory faers_counts does:
#   data <- faers_get(.object, field)
#   groups <- c("primaryid", .events)
#   data <- unique(data, by = groups)        # SELECT DISTINCT ...
#   data[, list(N = .N), keyby = .events]    # GROUP BY ... COUNT(*)
#
# The db path issues an equivalent single SQL statement and only pulls the
# *aggregated* result (tiny) into memory.  `.na.rm` maps to a WHERE clause.
#
# Consistency rule (design doc §3.5): an arbitrary `.fn` preprocessor cannot be
# translated to SQL.  When `.fn` is non-NULL the field is materialized and the
# original in-memory logic runs, preserving correctness at the cost of some
# memory advantage.  Only `.fn = NULL` takes the SQL fast path.

counts_db <- function(.object, .events = "soc_name", .fn = NULL, ...,
                      .field = "reac", .na.rm = FALSE) {
    assert_bool(.na.rm)

    if (!is.null(.fn)) {
        # Black-box preprocessor: fall back to the in-memory logic after
        # materializing the field (design doc §3.5).
        data <- faers_get(.object, field = .field)
        data <- rlang::as_function(.fn)(data, ...)
        if (!data.table::is.data.table(data)) {
            cli::cli_abort("{.fn .fn} must return a {.cls data.table}")
        }
        groups <- c("primaryid", .events)
        data <- unique(data, by = groups, cols = character())
        if (.na.rm) {
            keep <- !Reduce(`|`, lapply(data[, .SD, .SDcols = .events], is.na))
            data <- data[keep]
        }
        return(eval(substitute(
            data[, list(N = .N), keyby = .events],
            list(.events = .events)
        )))
    }

    # SQL fast path.  If the field's event columns come from the MedDRA
    # hierarchy (standardized reac/indi), join the small meddra lookup instead
    # of materializing the full field.
    meddra <- field_needs_meddra(.object, .field)
    if (meddra && !DB_MEDDRA_TABLE %in% DBI::dbListTables(.object@db@con)) {
        db_ingest_meddra(.object@db@con, .object)
    }
    q <- db_counts_query(.field, .events, na_rm = .na.rm, meddra = meddra)
    out <- DBI::dbGetQuery(.object@db@con, q)
    out <- data.table::as.data.table(out)
    # Byte-parity with the memory path (keyby + integer counts):
    #  - COUNT(DISTINCT) comes back BIGINT/double -> coerce to integer
    #  - keyby sets a data.table key on the event columns
    if (is.double(out$N)) data.table::set(out, j = "N", value = as.integer(out$N))
    data.table::setkeyv(out, .events)
    out
}
