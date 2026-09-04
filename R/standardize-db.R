# Database-mode standardization.
#
# Per design doc §3.6 (phase 1), standardize is treated as a materializing
# transform: for db-backed objects, each affected field (indi, reac) is pulled
# into memory one at a time, run through the *exact same* MedDRA mapping code
# as the memory path, then written back to DuckDB and released.  Because it
# reuses `meddra_standardize_pt` / `clean_indi_pt` / `clean_reac_pt`
# verbatim, the resulting identifiers are byte-for-byte identical to the
# memory path.  Only one field is resident at a time, so peak memory stays low
# even when MedDRA columns are added.

standardize_db <- function(object, meddra_path, add_smq = FALSE) {
    assert_string(meddra_path)
    meddra_data <- meddra(meddra_path, add_smq = add_smq, primary_soc = TRUE)

    meddra_cols <- c(
        "meddra_hierarchy_idx", "meddra_hierarchy_from",
        "meddra_code", "meddra_pt"
    )

    # ---- indi ----
    cli::cli_alert("standardize {.field Preferred Term} in indi")
    indi <- db_collect(object@db@con, "indi")
    indi <- data.table::as.data.table(indi)
    indi[, (meddra_cols) := meddra_standardize_pt(
        clean_indi_pt(indi_pt, meddra_data@hierarchy), # nolint
        meddra_data@hierarchy
    )]
    db_ingest(object@db@con, "indi", indi)
    rm(indi); gc()

    # ---- reac ----
    cli::cli_alert("standardize {.field Preferred Term} in reac")
    reac <- db_collect(object@db@con, "reac")
    reac <- data.table::as.data.table(reac)
    reac[, (meddra_cols) := meddra_standardize_pt(
        clean_reac_pt(pt, meddra_data@hierarchy), # nolint
        meddra_data@hierarchy
    )]
    db_ingest(object@db@con, "reac", reac)
    rm(reac); gc()

    object@meddra <- meddra_data
    object@standardization <- TRUE
    object
}
