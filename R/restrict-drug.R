# faers_reports_with: report-set constructor over drug exposures.
#
# A report-set builder: narrows a standardized FAERS object to reports whose
# drug records satisfy a combination of exposure conditions (multiple drug
# patterns, per-drug role constraints, AND/OR combination, exclusion of
# concomitant drugs). This is a set-construction tool for report counting,
# set operations and sub-population exploration; it is NOT by itself a
# drug-interaction / combination-toxicity signal analysis — scanning a subset
# built here lacks an unexposed comparator group (see the man page notes).
#
# Semantics (frozen):
# - Each pattern must hold on AT LEAST ONE drug row of the report; under
#   `.logic = "and"` two patterns MAY be satisfied by the SAME row (fixed-dose
#   combinations / multi-ingredient prod_ai entries).
# - `.roles` is a strictly-named list keyed by `.drugs` elements; drugs absent
#   from it are unrestricted. Unknown/duplicated names and invalid role codes
#   abort.
# - `.exclude_drugs` is report-level: any matching drug row (any role) removes
#   the report.
# - Matching semantics: grepl(pattern, x, ignore.case = TRUE). For the duckdb
#   backend a pattern is pushed down as regexp_matches(col, p, 'i') ONLY when
#   it is pushdown-safe (plain literal: letters/digits/space/underscore/hyphen);
#   otherwise the drug table is collected and matched in R, so both backends
#   agree byte-for-byte on every pattern.
# - Legal conditions with zero hits return an empty set (empty object /
#   character(0) / 0L) — this is a constructor, not a filter; an empty result
#   is a legitimate outcome.

# ---------- validation ----------

# Normalize + validate `.drugs`: non-empty character vector of regexes.
# Compiles each pattern so invalid regex aborts here.
.reports_with_check_drugs <- function(.drugs, call = rlang::caller_env()) {
    if (is.null(.drugs) || !is.character(.drugs) || !length(.drugs)) {
        cli::cli_abort(
            "{.arg .drugs} must be a non-empty character vector of drug name patterns",
            call = call
        )
    }
    .drugs <- trimws(.drugs)
    if (anyNA(.drugs) || any(!nzchar(.drugs))) {
        cli::cli_abort(
            "{.arg .drugs} must not contain NA or empty strings",
            call = call
        )
    }
    if (anyDuplicated(.drugs)) {
        cli::cli_abort(c(
            "{.arg .drugs} must not contain duplicated patterns",
            x = "Duplicated: {.val {unique(.drugs[duplicated(.drugs)])}}"
        ), call = call)
    }
    bad <- vapply(.drugs, function(p) {
        suppressWarnings(is.na(tryCatch(grepl(p, ""), error = function(e) NA)))
    }, logical(1L))
    if (any(bad)) {
        bad_val <- .drugs[which(bad)]
        cli::cli_abort(c(
            "Invalid regular expression in {.arg .drugs}",
            x = "{.val {bad_val}}"
        ), call = call)
    }
    .drugs
}

# Validate `.roles`: strictly-named list keyed by `.drugs` elements.
# Drugs absent from the list are unrestricted (no recycling).
.reports_with_check_roles <- function(.roles, .drugs, call = rlang::caller_env()) {
    if (is.null(.roles)) {
        return(list())
    }
    if (!is.list(.roles) || !length(.roles) ||
        is.null(names(.roles)) || any(!nzchar(names(.roles)))) {
        cli::cli_abort(
            "{.arg .roles} must be a fully named list like {.code list(warfarin = c(\"PS\", \"SS\"))}",
            call = call
        )
    }
    nms <- names(.roles)
    dup <- unique(nms[duplicated(nms)])
    if (length(dup)) {
        cli::cli_abort(c(
            "{.arg .roles} contains duplicated drug names",
            x = "{.val {dup}}"
        ), call = call)
    }
    unknown <- setdiff(nms, .drugs)
    if (length(unknown)) {
        cli::cli_abort(c(
            "{.arg .roles} refers to drug names not listed in {.arg .drugs}",
            x = "Unknown: {.val {unknown}}",
            i = "Role constraints may only be attached to drugs listed in {.arg .drugs}; unspecified drugs are unrestricted"
        ), call = call)
    }
    allowed <- c("PS", "SS", "C", "I")
    for (nm in nms) {
        roles <- .roles[[nm]]
        if (!is.character(roles) || !length(roles)) {
            cli::cli_abort(c(
                "{.arg .roles[[\"{nm}\"]]}} must be a non-empty character vector",
                i = "Valid role codes: {.val {allowed}}"
            ), call = call)
        }
        roles <- toupper(trimws(roles))
        bad <- setdiff(roles, allowed)
        if (length(bad)) {
            cli::cli_abort(c(
                "Invalid role code {.val {bad}} for drug {.val {nm}}",
                i = "Valid codes: {.val {allowed}} (primary suspect, secondary suspect, concomitant, interacting)"
            ), call = call)
        }
        .roles[[nm]] <- roles
    }
    .roles
}

# Validate `.exclude_drugs`: optional character vector of patterns.
.reports_with_check_exclude <- function(.exclude_drugs, call = rlang::caller_env()) {
    if (is.null(.exclude_drugs)) {
        return(character())
    }
    if (!is.character(.exclude_drugs) || !length(.exclude_drugs)) {
        cli::cli_abort(
            "{.arg .exclude_drugs} must be a character vector of drug patterns or NULL",
            call = call
        )
    }
    .exclude_drugs <- trimws(.exclude_drugs)
    if (anyNA(.exclude_drugs) || any(!nzchar(.exclude_drugs))) {
        cli::cli_abort(
            "{.arg .exclude_drugs} must not contain NA or empty strings",
            call = call
        )
    }
    unique(.exclude_drugs)
}

# Pushdown-safety: a pattern may go to duckdb regexp_matches(col, p, 'i')
# ONLY if it is a plain literal (letters/digits/space/underscore/hyphen).
# Anything else (metacharacters, classes, lookarounds) falls back to R-side
# grepl so both backends agree byte-for-byte (probed duckdb 1.5.5: RE2 bracket
# semantics differ from PCRE; lookarounds error out).
#
# An ASCII-only pattern is still NOT sufficient for identical case folding:
# RE2's 'i' flag applies simple Unicode folding (Kelvin sign U+212A "K" folds
# to "k") while R's TRE does not, so a non-ASCII drug name can match on
# duckdb but not in R. The SQL condition therefore also requires matched rows
# to be pure ASCII (`regexp_matches(col, '^[ -~]*$')`); the caller unions the
# remaining (non-ASCII) rows matched in R, where TRE folding defines the
# semantics for both backends.
.reports_with_pushdown_safe <- function(patterns) {
    all(grepl("^[A-Za-z0-9 _-]+$", patterns))
}

# ---------- shared matching core (memory backend) ----------

# Row-level fold matching with optional per-pattern role restriction. Returns
# a logical vector aligned with `drug`'s rows. R matching defines the
# semantics for the memory backend and the non-ASCII / regex SQL fallback.
.reports_with_row_hit <- function(drug, pattern, roles, drug_col) {
    out <- grepl(pattern, drug[[drug_col]], ignore.case = TRUE)
    if (length(roles)) out <- out & drug$role_cod %chin% roles
    out
}

# Build one report set per pattern before applying any report-level logic.
.reports_with_sets_mem <- function(drug, patterns, roles, drug_col) {
    lapply(patterns, function(p) {
        unique(drug$primaryid[
            .reports_with_row_hit(drug, p, roles[[p]], drug_col)
        ])
    })
}

# Every pattern's set already includes all matching rows, regardless of
# character encoding. Exclusions are also report sets from the complete
# drug table, and are applied only after AND/OR has been resolved.
.reports_with_combine <- function(id_sets, logic, exclude_sets) {
    ids <- if (identical(logic, "and")) Reduce(intersect, id_sets) else Reduce(union, id_sets)
    excluded <- Reduce(union, exclude_sets, init = character())
    setdiff(ids, excluded)
}

# ---------- db backend ----------

# SQL fragment: pattern restricted to its role codes. Literal patterns are
# pushed down (see .reports_with_pushdown_safe); everything else returns NULL
# and the caller falls back to the R path. Pushed-down matches are confined
# to pure-ASCII drug values so RE2 and TRE agree on the case folding; the
# caller re-matches every non-ASCII row in R and unions the results.
.reports_with_sql_cond <- function(pattern, roles, drug_col) {
    if (!.reports_with_pushdown_safe(pattern)) {
        return(NULL)
    }
    esc <- gsub("'", "''", pattern, fixed = TRUE)
    sql <- sprintf(
        "regexp_matches(%s, '%s', 'i') AND regexp_matches(%s, '^[ -~]*$')",
        db_qident(drug_col), esc, db_qident(drug_col)
    )
    if (length(roles)) {
        sql <- sprintf(
            "(%s AND role_cod IN (%s))",
            sql, paste0("'", roles, "'", collapse = ", ")
        )
    }
    sql
}

# Union SQL's ASCII hits and R's non-ASCII hits separately for EACH pattern.
# Combining the two halves after AND/exclusion would miss cross-half matches
# or reintroduce reports whose exclusion was in the other half.
.reports_with_sets_db <- function(con, patterns, roles, drug_col, nonascii) {
    nonascii_sets <- .reports_with_sets_mem(nonascii, patterns, roles, drug_col)
    lapply(seq_along(patterns), function(i) {
        cond <- .reports_with_sql_cond(patterns[[i]], roles[[patterns[[i]]]], drug_col)
        ascii_ids <- DBI::dbGetQuery(con, sprintf(
            "SELECT DISTINCT primaryid FROM %s WHERE %s",
            db_qident(db_field_table("drug")), cond
        ))$primaryid
        union(ascii_ids, nonascii_sets[[i]])
    })
}

# SQL fragment collecting the non-ASCII drug rows the pushed-down condition
# deliberately skipped: their folding is decided in R by TRE, which both
# backends treat as the baseline. Only (primaryid, drug value, role) rows that
# the ASCII guard excluded are fetched, so this stays cheap in practice.
.reports_with_nonascii_db <- function(con, drug_col) {
    DBI::dbGetQuery(
        con,
        sprintf(
            "SELECT primaryid, %s AS drug_value, role_cod FROM %s WHERE %s IS NOT NULL AND NOT regexp_matches(%s, '^[ -~]*$')",
            db_qident(drug_col), db_qident(db_field_table("drug")),
            db_qident(drug_col), db_qident(drug_col)
        )
    )
}

# ---------- main entry ----------

#' Build a report set from drug exposure conditions
#'
#' `faers_reports_with()` constructs a report sub-population from a
#' standardized FAERS object: reports whose drug records match a combination
#' of drug name patterns, optionally restricted to specific
#' [role codes][faers_phv_scan] and optionally excluding reports containing
#' other (concomitant) drugs.
#'
#' @details
#' ## What this function is — and is not
#'
#' `faers_reports_with()` is a **report set constructor**: it computes the set
#' of reports matching an exposure definition, for report counting, set
#' operations (AND/OR/exclusion) and sub-population exploration.
#'
#' It is **not** a drug-interaction or combination-toxicity analysis. Running
#' a disproportionality scan inside the returned subset uses the same exposure
#' definition for every drug, so there is no unexposed comparator group and
#' the 2x2 contingency table is not defined there. Keep set construction and
#' signal scanning as separate steps.
#'
#' ## Matching semantics
#'
#' Each pattern is a case-insensitive regular expression evaluated against
#' every drug row of a report (`grepl(..., ignore.case = TRUE)` semantics).
#' A report matches a pattern when **at least one** of its drug rows satisfies
#' it; under `.logic = "and"` all patterns must be satisfiable, and two
#' patterns **may be satisfied by the same drug row** (relevant for
#' fixed-dose combination products or `prod_ai` entries listing several
#' ingredients in one row).
#'
#' `.roles` restricts *which drug rows* may satisfy a pattern (e.g. only
#' primary-suspect rows). Unlisted drugs are unrestricted; role-constrained
#' and unrestricted rows are pooled per pattern (union within a drug).
#'
#' `.exclude_drugs` is report-level and role-agnostic: any drug row matching
#' an exclusion pattern removes the whole report, regardless of its role code.
#'
#' ## DuckDB pushdown
#'
#' For database-backed objects a pattern is pushed into SQL
#' (`regexp_matches(col, pattern, 'i')`) only when it is a plain literal
#' (letters, digits, spaces, underscores, hyphens). Any other pattern —
#' character classes, quantifiers, lookarounds — is matched in R after
#' collecting the drug table, so both backends always agree on the matched
#' report set. SQL and R matches are combined per pattern before applying
#' AND/OR logic and report-level exclusions. Report identifiers are collected
#' to combine these sets, including for `"count"` results.
#'
#' ## Empty results
#'
#' Legal conditions that match zero reports are not an error: an empty object,
#' `character(0)` or `0L` is returned. `faers_reports_with()` is a
#' constructor — an empty set is a legitimate result.
#'
#' @param .object A [FAERSascii] object. Standardization is required;
#'   de-duplication is recommended (duplicate drug rows do not double-count
#'   reports).
#' @param .drugs Character vector of drug name patterns (regular
#'   expressions, case-insensitive). Each must match at least one drug row of
#'   a report.
#' @param .roles A named list restricting the acceptable `role_cod` values
#'   per drug, e.g. `list(warfarin = c("PS", "SS"))`. Names must be elements
#'   of `.drugs`; drugs without an entry are unrestricted. Valid codes:
#'   `"PS"` (primary suspect), `"SS"` (secondary suspect), `"C"`
#'   (concomitant), `"I"` (interacting).
#' @param .logic `"and"` (default): every pattern must match. `"or"`: at
#'   least one pattern must match.
#' @param .exclude_drugs Optional character vector of patterns; reports with
#'   any drug row matching one of them are removed (any role).
#' @param .drug_col Optional name of the drug column to match against
#'   (default: `drugname`, falling back to `prod_ai`).
#' @param .return What to return: `"object"` (default, the filtered
#'   [FAERSascii]), `"ids"` (sorted unique `primaryid` character vector) or
#'   `"count"` (integer scalar).
#'
#' @return
#' * `.return = "object"`: the filtered [FAERSascii] object.
#' * `.return = "ids"`: sorted, unique `primaryid` character vector.
#' * `.return = "count"`: `integer(1)`, the number of matched reports.
#'
#' @examples
#' data <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
#' data <- faers_dedup(data)
#'
#' # reports containing Humalog in any role
#' n <- faers_reports_with(data, "humalog", .return = "count")
#'
#' # AND combination + per-drug role constraints + exclusion of a comedication
#' ids <- suppressMessages(faers_reports_with(
#'     data,
#'     c("humalog", "glyburide"),
#'     .roles = list(humalog = c("PS", "SS")),
#'     .exclude_drugs = "avandia",
#'     .return = "ids"
#' ))
#' length(ids)
#'
#' @seealso [faers_restrict()] for report-level (demographic) restrictions,
#'   [faers_phv_scan()] for signal detection.
#' @export
faers_reports_with <- function(.object, .drugs, .roles = NULL,
                               .logic = c("and", "or"),
                               .exclude_drugs = NULL, .drug_col = NULL,
                               .return = c("object", "ids", "count")) {
    .logic <- match.arg(.logic)
    .return <- match.arg(.return)
    .drugs <- .reports_with_check_drugs(.drugs)
    .roles <- .reports_with_check_roles(.roles, .drugs)
    .exclude_drugs <- .reports_with_check_exclude(.exclude_drugs)
    assert_s3_class(.object, "FAERSascii")
    if (!.object@standardization) {
        cli::cli_abort("{.arg .object} must be standardized using {.fn faers_standardize}")
    }

    drug_col <- .reports_with_pick_col(.object, .drug_col)

    # Match each inclusion and exclusion against the full logical drug
    # table before combining report sets. Only non-ASCII rows need collecting
    # on the literal SQL path; complex regexes use the existing R fallback.
    if (!is.null(.object@db) &&
        .reports_with_pushdown_safe(c(.drugs, .exclude_drugs))) {
        con <- .object@db@con
        nonascii <- data.table::as.data.table(.reports_with_nonascii_db(con, drug_col))
        data.table::setnames(nonascii, "drug_value", drug_col)
        id_sets <- .reports_with_sets_db(con, .drugs, .roles, drug_col, nonascii)
        exclude_sets <- .reports_with_sets_db(
            con, .exclude_drugs, list(), drug_col, nonascii
        )
    } else {
        drug <- faers_get(.object, "drug")
        id_sets <- .reports_with_sets_mem(drug, .drugs, .roles, drug_col)
        exclude_sets <- .reports_with_sets_mem(drug, .exclude_drugs, list(), drug_col)
    }
    ids <- .reports_with_combine(id_sets, .logic, exclude_sets)
    ids <- sort(unique(as.character(ids)))

    # Count the same per-pattern report sets used above, including non-ASCII
    # hits and deduplicating reports that match in both halves.
    per_hit <- vapply(id_sets, length, integer(1L))
    if (length(ids) || any(per_hit > 0L)) {
        cli::cli_inform(
            "{.fn faers_reports_with} matched {length(ids)} report{?s}: {.val {paste0(.drugs, \": \", per_hit)}}"
        )
    }

    # ---- return forms ----
    if (identical(.return, "count")) {
        return(as.integer(length(ids)))
    }
    if (identical(.return, "ids")) {
        return(as.character(ids))
    }
    # "object": db objects get a private deep copy BEFORE faers_keep, because
    # faers_keep_db() physically DELETEs rows from the shared connection.
    if (!is.null(.object@db)) {
        .object <- restrict_private_db_copy(.object)
    }
    faers_keep(.object, ids, invert = FALSE)
}

# Pick the drug name column: `drugname` when present, else `prod_ai`,
# unless the user names one explicitly. Errors when nothing usable exists
# (e.g. a quarter lacking prod_ai) instead of silently matching nothing.
.reports_with_pick_col <- function(.object, .drug_col = NULL) {
    if (!is.null(.object@db)) {
        # DESCRIBE is metadata-only: no table materialization on the db path.
        drug_cols <- DBI::dbGetQuery(
            .object@db@con, sprintf("DESCRIBE %s", db_field_table("drug"))
        )$column_name
    } else {
        drug_cols <- names(faers_get(.object, "drug"))
    }
    drug_col_nm <- .drug_col
    if (!is.null(.drug_col)) {
        if (!.drug_col %chin% drug_cols) {
            cli::cli_abort(c(
                "{.arg .drug_col} {.val {drug_col_nm}} not found in drug data",
                i = "Available: {.val {drug_cols}}"
            ))
        }
        return(.drug_col)
    }
    if ("drugname" %chin% drug_cols) {
        return("drugname")
    }
    if ("prod_ai" %chin% drug_cols) {
        cli::cli_inform(
            "{.field drugname} not found; falling back to {.field prod_ai} for drug matching"
        )
        return("prod_ai")
    }
    cli::cli_abort(c(
        "Neither {.val {\"drugname\"}} nor {.val {\"prod_ai\"}} found in drug data",
        i = "Pass {.arg .drug_col} explicitly to choose a column"
    ))
}