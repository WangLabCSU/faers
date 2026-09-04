# FAERSdb / FAERSdbTbl: database-backed storage for FAERS quarterly data.
#
# When a FAERSascii object is created with `database = "duckdb"`, its data is
# held on disk in a DuckDB database instead of in memory as data.table objects.
# The `@db` slot on the FAERS class carries the store (open connection +
# physical table names); `@data` holds one lightweight FAERSdbTbl proxy per
# field so the existing accessor contract (`faers_get`, `faers_mget`, ...)
# keeps working unchanged.  Only the truly needed fields / aggregations are
# ever pulled back into memory.
#
# Design note: this is the "option C (hybrid)" line from the integration
# design doc.  The default `database = "memory"` path never touches this code
# and behaves byte-for-byte identically to the original package.

# A lightweight proxy for a single field held in a DuckDB connection.
# `@con` points at the same (open) connection owned by the enclosing FAERSdb.
methods::setClass(
    "FAERSdbTbl",
    slots = list(
        con   = "ANY",      # the shared DuckDBConnection
        table = "character" # physical table name
    )
)

# The database-backed store attached to a FAERSascii object via `@db`.
# Holds one *open* connection for the object's lifetime; all per-field access
# reuses it, so fields are only materialized on demand.
methods::setClass(
    "FAERSdb",
    slots = list(
        con     = "ANY",       # an open DuckDBConnection
        path    = "character", # dbdir (file path, or ":memory:")
        version = "character", # duckdb package version used to create it
        tables  = "character"  # named char vector: field -> physical table name
    )
)

methods::setClassUnion("FAERSdbOrNULL", c("NULL", "FAERSdb"))

utils::globalVariables(c("duckdb", "DBI"))
