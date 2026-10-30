# faers 1.5.6

* Added `faers_reports_with()`, a report set constructor for drug exposure
  definitions: reports whose drug records match a combination of drug name
  patterns (case-insensitive regular expressions), optionally restricted per
  drug to specific role codes (`.roles`) and optionally excluding reports
  containing other drugs (`.exclude_drugs`, report-level and role-agnostic).
  `.logic = "and"/"or"` combines patterns; two patterns may be satisfied by
  the same drug row (fixed-dose combination products). Returns the filtered
  object, sorted `primaryid` vector or a count. On duckdb-backed objects
  plain-literal patterns are pushed down to SQL (`regexp_matches`), other
  patterns fall back to identical R matching. Per-drug matches are combined
  across ASCII and non-ASCII rows before AND/OR and report-level exclusions;
  object-mode returns build on a private copy, leaving the source database
  untouched. Per-pattern hit counts are reported as a message; empty result
  sets are a legitimate outcome, not an error.

* Added `faers_restrict()`, a report-level restriction helper that narrows a
  standardized FAERS object to a sub-population by reporting period
  (`.period`), date window (`.date_range` on `.date_field`, with partial
  dates like year-only `1999` excluded and counted by default, or kept with
  `.date_partial = "keep"`), sex, age interval (on the standardized
  `age_in_years`), age group and occurrence country. All enabled conditions
  are intersected; excluded reports (including missing or unparseable values)
  are counted per condition and reported in a single aggregated warning. The
  implementation converges onto `primaryid` and reuses `faers_keep()`, so
  both the memory and duckdb backends work and the
  `@standardization` / `@deduplication` flags are preserved.

* Added `.role_cod` to `faers_phv_scan()` and `faers_phv_scan_stratified()`:
  restrict the exposure definition to specific drug role codes (`PS` primary
  suspect, `SS` secondary suspect, `C` concomitant, `I` interacting) before
  counting. `NULL` (default) keeps every role, matching previous behaviour;
  the background `n` / `n_stratum` is unchanged either way. Implemented on
  both backends (memory pre-filter and SQL `role_cod IN (...)` pushdown).

* Added `faers_phv_scan()`, a whole-database signal scanner that enumerates
  all drug x event pairs, builds a contingency table for every pair and runs
  the requested disproportionality methods in one vectorized pass. It ships
  with an optimized backend for `database = "duckdb"` (single out-of-core SQL
  aggregation) and a chunked implementation for the memory backend, and
  supports optional drug whitelisting via `.drug_pattern`.

* Fixed `faers_phv_composite()`: hierarchy event types such as `soc_name`
  no longer abort with "object not found" (event columns are now resolved
  via `faers_get()`, which attaches the MedDRA hierarchy on both the memory
  and the duckdb backends), any standardized `reac` event column is accepted,
  an unused full-database `faers_counts()` call that crashed on duckdb-backed
  objects was removed, de-duplication is now enforced, and the undeclared
  `assertthat` usage was replaced with the internal assertion helpers.

* Added `faers_phv_scan_stratified()`, a stratified variant of the
  whole-database scanner that splits the database into strata defined by
  categorical `demo` columns (`sex`, `age_grp`, `occp_cod`, ...) and builds
  one 2x2 contingency table per `stratum x drug x event` combination against
  the stratum-internal background (`n_stratum`). Reports with missing strata
  values are grouped into an explicit `"Missing"` stratum by default
  (`.na_stratum = "keep"`) or dropped (`"drop"`), with a warning when a
  stratification column is missing for more than 30% of the reports. Both
  backends produce identical results, with the duckdb backend pushing the
  whole stratified aggregation into a single SQL query.

# faers 1.5.5

* Documented the DuckDB backend: `faers()`, `faers_parse()`, `faers_standardize()`
  and the downstream pipeline now accept `database = "duckdb"`, keeping the
  quarterly data out of memory in an on-disk DuckDB database with an identical
  analysis API (`faers_get`, `faers_counts`, `faers_phv_table`, ...).

# faers 1.5.1

* Fixed error in `handle_setopt(h, ...)` caused by unsupported option `multi_timeout`.
  Requires `curl` >= 6.0.0.

# faers 1.1.6

* fda_drugs() now directly use a fixed url to download the data

* faers_meta(internal = TRUE) will always use the cache data in the package. 

* Rename "gndr_cod" into "gender" for periods before 2014q2, and rename "sex" into "gender" for periods after or equal to 2014q2. 

* "sex" was added, which recoded any values other than "F" or "M" as `NA`.

# faers 1.1.4

* fix error when download failed

* meddra augment additional argument `primary_soc` to help save the full meddra data

# faers 0.99.0

* Initial Bioconductor submission.
