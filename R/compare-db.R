# DuckDB helpers for multi-group comparison (`faers_compare()`).
#
# Why this file exists: per-group event counts need the reports of a group as
# a set, but the public `faers_keep()` cannot be used here. On a duckdb-backed
# object `faers_keep_db()` physically DELETEs rows from the connection shared
# by every object built from the same database (see R/methods-FAERS.R), and the
# only safe alternative, `restrict_private_db_copy()` (R/restrict.R), deep
# copies the whole database once per group. A multi-group comparison would
# therefore pay one full copy per group, which is not acceptable for a
# multi-quarter on-disk database.
#
# Instead the group is registered as a TEMP table and the aggregation happens
# in a single read-only SELECT. The user's tables are never modified.
#
# The temp table name is generated per call rather than fixed. A fixed name is
# not safe: `dbWriteTable(temporary = TRUE)` does not shadow a same-named
# permanent table, it REPLACES it, so a fixed name would destroy whatever the
# user had called `__faers_compare_ids__` - and on an on-disk database the loss
# survives the session. The generated name is checked against the tables the
# connection already knows, so a collision cannot happen; the table is dropped
# again once the counts have been read.
#
# Results are aligned byte-for-byte with the in-memory path exactly like
# `counts_db()` (R/counts-db.R) does: BIGINT -> integer coercion, `setkeyv()`
# on the event columns. Callback preprocessing is shared by both backends in
# R/compare.R.

# Generates a name no table of this connection currently uses. Deterministic
# (no RNG, so the caller's random stream is untouched) and cheap: one
# `dbListTables()` per group next to the aggregate query that follows.
.compare_id_seq <- new.env(parent = emptyenv())
.compare_id_seq$n <- 0L

compare_id_table <- function(con) {
    existing <- DBI::dbListTables(con)
    repeat {
        .compare_id_seq$n <- .compare_id_seq$n + 1L
        nm <- sprintf(
            "__faers_compare_ids_%d_%d__", Sys.getpid(), .compare_id_seq$n
        )
        if (!nm %in% existing) {
            return(nm)
        }
    }
}

# Register `ids` as a temporary group table and return its name. A zero-length
# set is legal and yields a zero-row table, so empty groups need no special
# case downstream.
#
# `overwrite = FALSE` is deliberate rather than redundant: `compare_id_table()`
# already guarantees the name is free, and refusing to write over anything is
# the one behaviour that can never lose a user table. If the guarantee is ever
# broken the call errors out loudly instead of destroying data.
compare_ingest_ids <- function(con, ids) {
    tbl <- compare_id_table(con)
    DBI::dbWriteTable(
        con, tbl,
        data.frame(primaryid = as.character(ids)),
        overwrite = FALSE, temporary = TRUE
    )
    tbl
}

# Drop a temp table created by `compare_ingest_ids()`. Only ever called on a
# name this session generated, so it can never touch a user's table.
compare_drop_ids <- function(con, tbl) {
    DBI::dbExecute(con, sprintf("DROP TABLE IF EXISTS %s", db_qident(tbl)))
    invisible(con)
}

# One aggregate per event combination, restricted to the registered group.
# Mirrors `db_counts_query()` but scopes `COUNT(DISTINCT primaryid)` to the
# join against `id_table`. Every column reference is qualified because both
# sides of the join carry `primaryid`.
compare_counts_query <- function(field, events, id_table, na_rm = FALSE,
                                 meddra = FALSE, hierarchy = character()) {
    if (isTRUE(meddra)) {
        from <- sprintf(
            "%s AS f LEFT JOIN %s AS m ON f.meddra_hierarchy_idx = m.idx",
            db_field_table(field), DB_MEDDRA_TABLE
        )
        ev <- vapply(events, function(e) {
            if (e %chin% hierarchy) paste0("m.", db_qident(e)) else paste0("f.", db_qident(e))
        }, character(1L))
    } else {
        from <- sprintf("%s AS f", db_field_table(field))
        ev <- vapply(events, function(e) paste0("f.", db_qident(e)), character(1L))
    }
    join_ids <- sprintf(
        "INNER JOIN %s AS g ON f.primaryid = g.primaryid", db_qident(id_table)
    )
    where <- ""
    if (isTRUE(na_rm)) {
        where <- paste("WHERE", paste(sprintf("%s IS NOT NULL", ev), collapse = " AND "))
    }
    select <- if (length(ev) == 1L) {
        event_name <- if (nzchar(events[[1L]])) db_qident(events[[1L]]) else "event"
        sprintf("SELECT %s AS %s, COUNT(DISTINCT f.primaryid) AS N", ev, event_name)
    } else {
        sprintf("SELECT %s, COUNT(DISTINCT f.primaryid) AS N", paste(ev, collapse = ", "))
    }
    sprintf(
        "%s FROM %s %s %s GROUP BY %s",
        select, from, join_ids, where, paste(ev, collapse = ", ")
    )
}

# Database-mode per-group counting. Used by `faers_compare()` for every group
# on a duckdb-backed object; the in-memory backend keeps using
# `faers_keep()` + `faers_counts()`.
compare_counts_db <- function(.object, .ids, .events,
                              .field = "reac", .na.rm = FALSE) {
    assert_bool(.na.rm)
    con <- .object@db@con

    meddra <- field_needs_meddra(.object, .field)
    if (meddra && !DB_MEDDRA_TABLE %in% DBI::dbListTables(con)) {
        db_ingest_meddra(con, .object)
    }
    tbl <- compare_ingest_ids(con, .ids)
    on.exit(compare_drop_ids(con, tbl), add = TRUE)
    out <- data.table::as.data.table(DBI::dbGetQuery(
        con,
        compare_counts_query(
            .field, .events, id_table = tbl,
            na_rm = .na.rm, meddra = meddra,
            hierarchy = db_meddra_columns(.object)
        )
    ))
    # Byte-parity with the memory path (see counts_db()).
    if (is.double(out$N)) data.table::set(out, j = "N", value = as.integer(out$N))
    data.table::setkeyv(out, .events)
    out
}

# Narrow projection of the `drug` field used by the single-row combination
# diagnosis: `primaryid`, every drug column a group matches against, and
# `role_cod`. On duckdb only these columns are pulled instead of the whole
# drug table.
compare_drug_rows <- function(.object, cols) {
    cols <- setdiff(unique(cols), c("primaryid", "role_cod"))
    if (is.null(.object@db)) {
        drug <- faers_get(.object, "drug")
        return(drug[, c("primaryid", cols, "role_cod"), with = FALSE])
    }
    data.table::as.data.table(DBI::dbGetQuery(
        .object@db@con,
        sprintf(
            "SELECT %s FROM %s",
            paste(
                c("primaryid", vapply(cols, db_qident, character(1L)), "role_cod"),
                collapse = ", "
            ),
            db_field_table("drug")
        )
    ))
}
