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
# in a single read-only SELECT. The user's tables are never modified: only
# `__faers_compare_ids__` is created (and overwritten) in the session.
#
# Results are aligned byte-for-byte with the in-memory path exactly like
# `counts_db()` (R/counts-db.R) does: BIGINT -> integer coercion, `setkeyv()`
# on the event columns, and a materialize-and-fall-back path for an arbitrary
# `.fn` preprocessor that cannot be translated to SQL.

.COMPARE_ID_TABLE <- "__faers_compare_ids__"

# Register `ids` as the temporary group table. A zero-length set is legal and
# yields a zero-row table, so empty groups need no special case downstream.
compare_ingest_ids <- function(con, ids) {
    DBI::dbWriteTable(
        con, .COMPARE_ID_TABLE,
        data.frame(primaryid = as.character(ids)),
        overwrite = TRUE, temporary = TRUE
    )
    invisible(con)
}

# One aggregate per event combination, restricted to the registered group.
# Mirrors `db_counts_query()` but scopes `COUNT(DISTINCT primaryid)` to the
# join against `__faers_compare_ids__`. Every column reference is qualified
# because both sides of the join carry `primaryid`.
compare_counts_query <- function(field, events, na_rm = FALSE, meddra = FALSE,
                                 hierarchy = character()) {
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
        "INNER JOIN %s AS g ON f.primaryid = g.primaryid", .COMPARE_ID_TABLE
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

# Collect the rows of `field` that belong to `ids`, with the MedDRA hierarchy
# attached for `reac`/`indi` so the result matches `faers_get()` exactly. Used
# only for the black-box `.fn` fallback.
compare_collect_field <- function(.object, .ids, .field) {
    con <- .object@db@con
    compare_ingest_ids(con, .ids)
    out <- data.table::as.data.table(DBI::dbGetQuery(
        con,
        sprintf(
            "SELECT f.* FROM %s AS f INNER JOIN %s AS g ON f.primaryid = g.primaryid",
            db_field_table(.field), .COMPARE_ID_TABLE
        )
    ))
    if (.object@standardization && any(.field == c("indi", "reac"))) {
        out <- faers_add_meddra(out, .object@meddra@hierarchy)
    }
    out
}

# Database-mode per-group counting. Used by `faers_compare()` for every group
# on a duckdb-backed object; the in-memory backend keeps using
# `faers_keep()` + `faers_counts()`.
compare_counts_db <- function(.object, .ids, .events, .fn = NULL, ...,
                              .field = "reac", .na.rm = FALSE) {
    assert_bool(.na.rm)
    con <- .object@db@con

    if (!is.null(.fn)) {
        # Black-box preprocessor: same rule as `counts_db()` (design doc §3.5) -
        # materialize the field for the group and run the in-memory logic.
        data <- compare_collect_field(.object, .ids, .field)
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

    meddra <- field_needs_meddra(.object, .field)
    if (meddra && !DB_MEDDRA_TABLE %in% DBI::dbListTables(con)) {
        db_ingest_meddra(con, .object)
    }
    compare_ingest_ids(con, .ids)
    out <- data.table::as.data.table(DBI::dbGetQuery(
        con,
        compare_counts_query(
            .field, .events,
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
