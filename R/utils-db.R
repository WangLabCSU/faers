# Low-level DuckDB helpers used by the database-backed FAERS path.
#
# These are intentionally small and dependency-light: they talk to the DuckDB
# connection via DBI and return plain R types.  The S4 layer (FAERSdb) and the
# pipeline functions (parse / counts / phv / materialize) build on top.

# Is DuckDB available in this session? All db-mode functions must guard on this
# before touching duckdb so that the default (memory) path never requires it.
db_available <- function() {
    requireNamespace("duckdb", quietly = TRUE) && requireNamespace("DBI", quietly = TRUE)
}

# Assert that DuckDB is installed, with a helpful message otherwise.
db_require <- function() {
    if (!db_available()) {
        cli::cli_abort(c(
            "The {.emph duckdb} backend is not available on this machine.",
            i = "Install it with {.code install.packages(\"duckdb\")} to use {.arg database = \"duckdb\"}.",
            i = "The default {.arg database = \"memory\"} path does not need it."
        ))
    }
    invisible(TRUE)
}

# Normalize the db path argument.  A file path (default) or ":memory:".
db_normalize_path <- function(db_path = NULL) {
    if (is.null(db_path) || identical(db_path, ":memory:")) {
        return(":memory:")
    }
    assert_string(db_path, allow_empty = FALSE)
    db_path
}

# Open a DuckDB connection to `db_path`. Returns a DuckDBConnection.
db_connect <- function(db_path = NULL) {
    db_require()
    path <- db_normalize_path(db_path)
    duckdb::dbConnect(duckdb::duckdb(), dbdir = path)
}

# Gracefully disconnect a connection (shuts the instance down to release the
# file handle).
db_disconnect <- function(con) {
    if (!is.null(con)) {
        try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    }
    invisible(NULL)
}

# Deterministic physical table name for a FAERS field.
db_field_table <- function(field) {
    paste0("faers_", field)
}

# The full named mapping of field -> physical table name used throughout.
db_field_tables <- function() {
    stats::setNames(
        vapply(FAERS_ASCII_FILE_FIELDS, db_field_table, character(1L)),
        FAERS_ASCII_FILE_FIELDS
    )
}

# Generate a fresh on-disk database path (combine/materialize make their own
# DB so the source connections are untouched).
db_new_path <- function() {
    file.path(tempdir(), sprintf("faers_%s.duckdb", digest_short(stats::runif(1L))))
}

# A short stable hash used to namespace temp objects / file names without
# pulling in a heavy hashing dependency.
digest_short <- function(x) {
    h <- as.character(utils::packageVersion("base"))
    # fall back to a simple hash so we never depend on extra packages
    s <- paste(x, collapse = "")
    n <- 0
    for (ch in utf8ToInt(s)) n <- (n * 31L + ch) %% 1000000007L
    as.character(abs(n))
}

# Write an in-memory data.table into the DuckDB `faers_<field>` table.
db_ingest <- function(con, field, dt) {
    table_name <- db_field_table(field)
    DBI::dbWriteTable(
        con, name = table_name, value = as.data.frame(dt),
        overwrite = TRUE, temporary = FALSE
    )
    invisible(table_name)
}

# Collect a field's table into a data.table.
db_collect <- function(con, field) {
    db_collect_table(con, db_field_table(field))
}

# Collect an already-prefixed physical table into a data.table.
db_collect_table <- function(con, table_name) {
    out <- DBI::dbGetQuery(con, sprintf("SELECT * FROM %s", table_name))
    data.table::as.data.table(out)
}

# Count rows in a field's table (no materialization).
db_nrow <- function(con, field) {
    table_name <- db_field_table(field)
    DBI::dbGetQuery(
        con, sprintf("SELECT COUNT(*) AS n FROM %s", table_name)
    )$n[[1L]]
}

# DuckDB SQL-quote an identifier (table/column name).
db_qident <- function(x) {
    DBI::dbQuoteIdentifier(DBI::ANSI(), x)
}

# Materialize every field of a db-backed FAERSascii into in-memory
# data.tables, returning a standard (memory) FAERSascii with `@db` cleared.
# Used for backwards-compatible inspection (`faers_data`, `$`, `[[`, `[`) and
# at the end of a pipeline when the caller wants the plain object.
materialize <- function(object) {
    if (is.null(object@db)) {
        return(object)
    }
    con <- object@db@con
    fields <- names(object@db@tables)
    data <- lapply(fields, function(field) db_collect(con, field))
    methods::new(
        "FAERSascii",
        data = stats::setNames(data, fields),
        deletedCases = methods::slot(object, "deletedCases"),
        year = object@year, quarter = object@quarter,
        standardization = object@standardization,
        deduplication = object@deduplication,
        format = object@format,
        db = NULL
    )
}

# Name of the small MedDRA hierarchy lookup table (NOT a FAERS file field).
DB_MEDDRA_TABLE <- "faers_meddra"

# Ingest the MedDRA hierarchy lookup into the connection as a small table
# (87k rows) indexed by `idx`.  The in-memory path does
# `data[, names(hierarchy) := hierarchy[.__idx__.]]` — a pure vectorized lookup
# by `meddra_hierarchy_idx` — so this table lets the SQL fast path reproduce
# that exactly for standardized reac/indi without materializing the full field.
db_ingest_meddra <- function(con, object) {
    hier <- object@meddra@hierarchy
    idxs <- seq_along(hier[[1L]])
    tbl <- data.frame(
        idx = idxs,
        lapply(hier, function(v) v[idxs]),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    DBI::dbWriteTable(
        con, name = DB_MEDDRA_TABLE, value = tbl,
        overwrite = TRUE, temporary = FALSE
    )
    invisible(DB_MEDDRA_TABLE)
}

# Does counting/logic on this field need the MedDRA hierarchy join?  Only for
# standardized reac/indi, where faers_get enriches the field with the hierarchy
# columns (soc_name, pt_name, ...) before counting.
field_needs_meddra <- function(object, field) {
    isTRUE(object@standardization) && any(field == c("reac", "indi"))
}

# Build a SQL string that counts DISTINCT primaryid per event group, with an
# optional WHERE clause (used for .na.rm and for event-set filtering).
#
# When `meddra = TRUE` the field table is joined to the small MedDRA lookup so
# that hierarchy-derived event columns (soc_name, pt_name, ...) can be grouped
# exactly as the in-memory path derives them — without materializing the field.
#
# The output column for a single event is named after the event itself,
# matching the memory path's keyby column name.
db_counts_query <- function(field, events, na_rm = FALSE, meddra = FALSE) {
    tbl <- db_field_table(field)
    join <- ""
    if (isTRUE(meddra)) {
        tbl <- sprintf("%s AS f LEFT JOIN %s AS m ON f.meddra_hierarchy_idx = m.idx",
            tbl, DB_MEDDRA_TABLE)
        ev <- vapply(events, function(e) paste0("m.", db_qident(e)), character(1L))
        primaryid_sql <- "f.primaryid"
    } else {
        ev <- vapply(events, db_qident, character(1L))
        primaryid_sql <- "primaryid"
    }
    where <- ""
    if (isTRUE(na_rm)) {
        conds <- sprintf("%s IS NOT NULL", ev)
        where <- paste("WHERE", paste(conds, collapse = " AND "))
    }
    if (length(ev) == 1L) {
        event_name <- if (nzchar(events[[1L]])) db_qident(events[[1L]]) else "event"
        sprintf(
            "SELECT %s AS %s, COUNT(DISTINCT %s) AS N FROM %s %s GROUP BY %s",
            ev, event_name, primaryid_sql, tbl, where, ev
        )
    } else {
        sprintf(
            "SELECT %s FROM %s %s GROUP BY %s",
            paste(c(ev, sprintf("COUNT(DISTINCT %s) AS N", primaryid_sql)),
                collapse = ", "),
            tbl, where, paste(ev, collapse = ", ")
        )
    }
}
