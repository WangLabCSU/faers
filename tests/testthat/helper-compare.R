compare_db_twin <- function(mem_obj) {
    ns <- asNamespace("faers")
    con <- ns$db_connect(":memory:")
    tables <- ns$db_field_tables()
    for (f in names(mem_obj@data)) ns$db_ingest(con, f, mem_obj@data[[f]])
    ns$db_ingest_meddra(con, mem_obj)
    proxies <- stats::setNames(
        lapply(names(mem_obj@data), function(f) {
            methods::new("FAERSdbTbl", con = con, table = tables[[f]])
        }),
        names(mem_obj@data)
    )
    methods::new("FAERSascii",
        data = proxies,
        deletedCases = methods::slot(mem_obj, "deletedCases"),
        year = mem_obj@year, quarter = mem_obj@quarter,
        standardization = mem_obj@standardization,
        deduplication = mem_obj@deduplication,
        format = mem_obj@format,
        meddra = mem_obj@meddra,
        db = methods::new("FAERSdb",
            con = con, path = ":memory:",
            version = as.character(packageVersion("duckdb")),
            tables = tables
        )
    )
}

compare_std_data <- function() {
    data <- readRDS(system.file("extdata", "standardized_data.rds", package = "faers"))
    suppressMessages(faers_dedup(data))
}

compare_quiet <- function(expr) {
    suppressWarnings(suppressMessages(expr))
}
