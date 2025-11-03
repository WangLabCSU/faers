#' Format Pharmacovigilance Signal Detection Results
#' 
#' This function processes raw pharmacovigilance signal detection results and
#' generates formatted tables suitable for reporting and analysis. It provides
#' flexible options for method selection, significance filtering, and output formatting.
#'
#' @param pv_data A data frame containing pharmacovigilance signal detection results.
#' @param term_col Character. Column name for medical terminology (e.g., "soc_name", "pt_name").
#'   Default is "soc_name".
#' @param methods Character vector. Signal detection methods to include. Options include:
#'   "ror", "prr", "bcpnn_norm", "bcpnn_mcmc", "oe_ratio", "odds_ratio", "ebgm".
#'   Default is NULL (include all available methods).
#' @param show_significant_only Logical. Show only statistically significant signals?
#'   Default is FALSE.
#' @param clear_format Logical. Use clear format with standardized column names?
#'   Default is FALSE (keep original column names).
#' @param significance_criteria List. Custom significance criteria for each method.
#'   If NULL, uses standard criteria. Default is NULL.
#' @param min_cases Integer. Minimum number of cases required for inclusion.
#'   Default is 0 (no minimum).
#' @param sort_by Character. Column name to sort results by. Default is "ror".
#' @param descending Logical. Sort in descending order? Default is TRUE.
#' @param max_rows Integer. Maximum number of rows to return. Default is NULL (all rows).
#' @param drop_na Logical. Remove rows with NA in term column? Default is TRUE.
#' @param remove_zero_cases Logical. Remove rows with zero cases? Default is TRUE.
#'
#' @return A data frame with formatted pharmacovigilance results. The structure
#'   depends on the clear_format parameter:
#'   - If clear_format = FALSE: Returns selected columns with original names
#'   - If clear_format = TRUE: Returns standardized columns with clear names
#'
#' @details
#' The function supports the following signal detection methods with standard significance criteria:
#' \itemize{
#'   \item \strong{ROR}: ROR > 2 & ROR_CI_low > 1
#'   \item \strong{PRR}: PRR ≥ 2 & χ² ≥ 4 & cases ≥ 3
#'   \item \strong{BCPNN (norm)}: IC > 0 & IC025 > 0
#'   \item \strong{BCPNN (mcmc)}: IC > 0 & IC025 > 0  
#'   \item \strong{OE Ratio}: OE > 2 & OE_CI_low > 1
#'   \item \strong{Odds Ratio}: OR > 2 & OR_CI_low > 1
#'   \item \strong{EBGM}: EBGM > 2 & EBGM05 > 1
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage - all methods, original column names
#' basic_result <- format_pv_results(Immunity_signals)
#'
#' # Remove rows where the case count is 0.
#' filtered_result <- format_pv_results(
#'   Immunity_signals,
#'   remove_zero_cases = TRUE
#' )
#'
#' # Comprehensive filtering: zero cases, minimum case count, and significant signals.
#' comprehensive_result <- format_pv_results(
#'   Immunity_signals,
#'   remove_zero_cases = TRUE,
#'   min_cases = 3,
#'   show_significant_only = TRUE
#' )
#' }
#'
#' @export
format_pv_results <- function(pv_data,
                              term_col = "soc_name",
                              methods = NULL,
                              show_significant_only = FALSE,
                              clear_format = FALSE,
                              significance_criteria = NULL,
                              min_cases = 0,
                              sort_by = "ror",
                              descending = TRUE,
                              max_rows = NULL,
                              drop_na = TRUE,
                              remove_zero_cases = TRUE) {
  
  # Step 1: Input validation and preprocessing
  validated_data <- validate_and_preprocess_inputs(
    pv_data, term_col, sort_by, drop_na, min_cases, remove_zero_cases
  )
  
  # Step 2: Method selection and configuration
  method_config <- configure_methods(validated_data, methods)
  
  # Step 3: Significance filtering (if requested)
  if (show_significant_only) {
    validated_data <- apply_significance_filter(
      validated_data, method_config, significance_criteria
    )
    if (nrow(validated_data) == 0) {
      cli::cli_alert_warning("No significant signals found with the specified criteria.")
      return(NULL)
    }
  }
  
  # Step 4: Data sorting and row limiting
  processed_data <- sort_and_limit_data(
    validated_data, sort_by, descending, max_rows
  )
  
  # Step 5: Format final output
  final_table <- create_final_output(
    processed_data, term_col, method_config, clear_format
  )
  
  # Add metadata attributes
  final_table <- add_metadata_attributes(
    final_table, term_col, method_config, show_significant_only, clear_format
  )
  
  return(final_table)
}

# Input Validation and Preprocessing ------------------------------------------

validate_and_preprocess_inputs <- function(pv_data, term_col, sort_by, drop_na, min_cases, remove_zero_cases) {
  
  # Check for required base columns
  required_cols <- c(term_col, "a")
  missing_cols <- setdiff(required_cols, names(pv_data))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Missing required columns: {paste(missing_cols, collapse = ', ')}")
  }
  
  # Verify sort column exists
  if (!sort_by %in% names(pv_data)) {
    cli::cli_abort("Sort column '{sort_by}' not found in data")
  }
  
  # Convert to data.frame for consistent processing
  if (inherits(pv_data, "data.table")) {
    pv_data <- as.data.frame(pv_data)
    cli::cli_alert_info("Converted data.table to data.frame for processing")
  }
  
  # Remove NA values in term column if requested
  if (drop_na) {
    initial_rows <- nrow(pv_data)
    pv_data <- pv_data[!is.na(pv_data[[term_col]]), ]
    removed_rows <- initial_rows - nrow(pv_data)
    if (removed_rows > 0) {
      cli::cli_alert_info("Removed {removed_rows} rows with NA in {term_col}")
    }
  }
  
  # NEW: Remove rows with zero cases if requested
  if (remove_zero_cases) {
    initial_rows <- nrow(pv_data)
    pv_data <- pv_data[pv_data$a > 0, ]
    zero_case_rows <- initial_rows - nrow(pv_data)
    if (zero_case_rows > 0) {
      cli::cli_alert_info("Removed {zero_case_rows} rows with zero cases")
    }
  }
  
  # Apply minimum cases filter (after zero case removal)
  if (min_cases > 0) {
    initial_rows <- nrow(pv_data)
    pv_data <- pv_data[pv_data$a >= min_cases, ]
    filtered_rows <- initial_rows - nrow(pv_data)
    if (filtered_rows > 0) {
      cli::cli_alert_info("Filtered {filtered_rows} rows with fewer than {min_cases} cases")
    }
  }
  
  cli::cli_alert_success("Input validation completed: {nrow(pv_data)} rows remaining")
  return(pv_data)
}

# Method Configuration --------------------------------------------------------

configure_methods <- function(data, selected_methods = NULL) {
  
  # Define all available methods and their required columns
  method_definitions <- list(
    ror = list(
      display_name = "ROR",
      required_cols = c("ror", "ror_ci_low", "ror_ci_high"),
      standard_criteria = function(df) df$ror > 2 & df$ror_ci_low > 1
    ),
    prr = list(
      display_name = "PRR", 
      required_cols = c("prr", "prr_ci_low", "prr_ci_high", "chisq", "chisq_pvalue"),
      standard_criteria = function(df) df$prr >= 2 & df$chisq >= 4 & df$a >= 3
    ),
    bcpnn_norm = list(
      display_name = "BCPNN (Norm)",
      required_cols = c("bcpnn_norm_ic", "bcpnn_norm_ic_ci_low", "bcpnn_norm_ic_ci_high"),
      standard_criteria = function(df) df$bcpnn_norm_ic > 0 & df$bcpnn_norm_ic_ci_low > 0
    ),
    bcpnn_mcmc = list(
      display_name = "BCPNN (MCMC)",
      required_cols = c("bcpnn_mcmc_ic", "bcpnn_mcmc_ic_ci_low", "bcpnn_mcmc_ic_ci_high"),
      standard_criteria = function(df) df$bcpnn_mcmc_ic > 0 & df$bcpnn_mcmc_ic_ci_low > 0
    ),
    oe_ratio = list(
      display_name = "OE Ratio",
      required_cols = c("oe_ratio", "oe_ratio_ci_low", "oe_ratio_ci_high"),
      standard_criteria = function(df) df$oe_ratio > 2 & df$oe_ratio_ci_low > 1
    ),
    odds_ratio = list(
      display_name = "Odds Ratio",
      required_cols = c("odds_ratio", "odds_ratio_ci_low", "odds_ratio_ci_high", "fisher_pvalue"),
      standard_criteria = function(df) df$odds_ratio > 2 & df$odds_ratio_ci_low > 1
    ),
    ebgm = list(
      display_name = "EBGM",
      required_cols = c("ebgm", "ebgm_ci_low", "ebgm_ci_high"),
      standard_criteria = function(df) df$ebgm > 2 & df$ebgm_ci_low > 1
    )
  )
  
  # Identify available methods based on column presence
  available_methods <- names(method_definitions)[
    sapply(method_definitions, function(method) {
      all(method$required_cols %in% names(data))
    })
  ]
  
  if (length(available_methods) == 0) {
    cli::cli_abort("No valid signal detection methods found in the data")
  }
  
  # Determine which methods to use
  if (is.null(selected_methods)) {
    used_methods <- available_methods
    cli::cli_alert_info("Using all available methods: {paste(available_methods, collapse = ', ')}")
  } else {
    # Validate selected methods
    invalid_methods <- setdiff(selected_methods, names(method_definitions))
    if (length(invalid_methods) > 0) {
      cli::cli_alert_warning("Invalid methods specified: {paste(invalid_methods, collapse = ', ')}")
      selected_methods <- setdiff(selected_methods, invalid_methods)
    }
    
    # Check availability of selected methods
    available_selected <- intersect(selected_methods, available_methods)
    unavailable_methods <- setdiff(selected_methods, available_methods)
    
    if (length(unavailable_methods) > 0) {
      cli::cli_alert_warning(
        "Methods not available (missing columns): {paste(unavailable_methods, collapse = ', ')}"
      )
    }
    
    if (length(available_selected) == 0) {
      cli::cli_alert_warning("No selected methods available. Using all available methods.")
      used_methods <- available_methods
    } else {
      used_methods <- available_selected
    }
  }
  
  cli::cli_alert_success("Configured methods: {paste(used_methods, collapse = ', ')}")
  
  return(list(
    definitions = method_definitions,
    used_methods = used_methods,
    available_methods = available_methods
  ))
}

# Significance Filtering ------------------------------------------------------

apply_significance_filter <- function(data, method_config, custom_criteria = NULL) {
  
  if (length(method_config$used_methods) == 0) {
    return(data)
  }
  
  # Create significance filters for each method
  significance_filters <- list()
  
  for (method in method_config$used_methods) {
    method_def <- method_config$definitions[[method]]
    
    # Use custom criteria if provided, otherwise standard criteria
    if (!is.null(custom_criteria) && method %in% names(custom_criteria)) {
      sig_test <- custom_criteria[[method]]
    } else {
      sig_test <- method_def$standard_criteria
    }
    
    significance_filters[[method]] <- sig_test(data)
  }
  
  # Combine filters: keep rows that are significant by at least one method
  if (length(significance_filters) > 0) {
    combined_filter <- Reduce(`|`, significance_filters)
    filtered_data <- data[combined_filter, ]
    
    # Count significance by method for reporting
    sig_counts <- sapply(significance_filters, sum)
    total_sig <- sum(combined_filter)
    
    cli::cli_alert_info(
      "Significance filter applied: {total_sig} signals significant by at least one method"
    )
    
    if (total_sig > 0) {
      sig_details <- paste(
        names(sig_counts), "(", sig_counts, ")", 
        sep = "", collapse = ", "
      )
      cli::cli_alert_info("Signals by method: {sig_details}")
    }
    
    return(filtered_data)
  }
  
  return(data)
}

# Data Processing -------------------------------------------------------------

sort_and_limit_data <- function(data, sort_by, descending, max_rows) {
  
  # Sort data
  if (descending) {
    sorted_data <- data[order(-data[[sort_by]]), ]
  } else {
    sorted_data <- data[order(data[[sort_by]]), ]
  }
  
  # Limit rows if specified
  if (!is.null(max_rows) && max_rows > 0 && max_rows < nrow(sorted_data)) {
    sorted_data <- head(sorted_data, max_rows)
    cli::cli_alert_info("Limited to {max_rows} rows")
  }
  
  return(sorted_data)
}

# Output Formatting -----------------------------------------------------------

create_final_output <- function(data, term_col, method_config, clear_format) {
  
  if (clear_format) {
    return(create_clear_format_table(data, term_col, method_config))
  } else {
    return(create_original_format_table(data, term_col, method_config))
  }
}

create_clear_format_table <- function(data, term_col, method_config) {
  
  term_label <- get_standardized_term_label(term_col)
  
  # Start with basic columns
  result <- data.frame(
    Terminology = data[[term_col]],
    Cases = data[["a"]],
    stringsAsFactors = FALSE
  )
  colnames(result)[1:2] <- c(term_label, "Cases")
  
  # Add selected method columns with clear names
  for (method in method_config$used_methods) {
    method_cols <- get_clear_format_columns(method, data)
    result <- cbind(result, method_cols)
  }
  
  return(result)
}

create_original_format_table <- function(data, term_col, method_config) {
  
  # Collect all columns to keep
  cols_to_keep <- c(term_col, "a")
  
  for (method in method_config$used_methods) {
    method_def <- method_config$definitions[[method]]
    cols_to_keep <- c(cols_to_keep, method_def$required_cols)
  }
  
  # Keep only existing columns and remove duplicates
  existing_cols <- unique(intersect(cols_to_keep, names(data)))
  
  if (length(existing_cols) == 0) {
    cli::cli_abort("No selected columns found in the data")
  }
  
  result <- data[, existing_cols, drop = FALSE]
  
  return(result)
}

# Helper Functions ------------------------------------------------------------

get_standardized_term_label <- function(term_col) {
  labels <- c(
    "soc_name" = "System Organ Class",
    "hlgt_name" = "High Level Group Term", 
    "hlt_name" = "High Level Term",
    "pt_name" = "Preferred Term",
    "llt_name" = "Lowest Level Term"
  )
  
  if (term_col %in% names(labels)) {
    return(labels[[term_col]])
  } else {
    return(tools::toTitleCase(gsub("_", " ", term_col)))
  }
}

get_clear_format_columns <- function(method, data) {
  
  column_definitions <- list(
    ror = data.frame(
      ROR = sprintf("%.2f", data$ror),
      ROR_CI = sprintf("(%.2f-%.2f)", data$ror_ci_low, data$ror_ci_high),
      stringsAsFactors = FALSE
    ),
    prr = data.frame(
      PRR = sprintf("%.2f", data$prr),
      PRR_CI = sprintf("(%.2f-%.2f)", data$prr_ci_low, data$prr_ci_high),
      PRR_Chi2 = sprintf("%.1f", data$chisq),
      stringsAsFactors = FALSE
    ),
    bcpnn_norm = data.frame(
      BCPNN_Norm = sprintf("%.2f", data$bcpnn_norm_ic),
      BCPNN_Norm_CI = sprintf("(%.2f-%.2f)", data$bcpnn_norm_ic_ci_low, data$bcpnn_norm_ic_ci_high),
      stringsAsFactors = FALSE
    ),
    bcpnn_mcmc = data.frame(
      BCPNN_MCMC = sprintf("%.2f", data$bcpnn_mcmc_ic),
      BCPNN_MCMC_CI = sprintf("(%.2f-%.2f)", data$bcpnn_mcmc_ic_ci_low, data$bcpnn_mcmc_ic_ci_high),
      stringsAsFactors = FALSE
    ),
    oe_ratio = data.frame(
      OE_Ratio = sprintf("%.2f", data$oe_ratio),
      OE_CI = sprintf("(%.2f-%.2f)", data$oe_ratio_ci_low, data$oe_ratio_ci_high),
      stringsAsFactors = FALSE
    ),
    odds_ratio = data.frame(
      Odds_Ratio = sprintf("%.2f", data$odds_ratio),
      OR_CI = sprintf("(%.2f-%.2f)", data$odds_ratio_ci_low, data$odds_ratio_ci_high),
      stringsAsFactors = FALSE
    ),
    ebgm = data.frame(
      EBGM = sprintf("%.2f", data$ebgm),
      EBGM_CI = sprintf("(%.2f-%.2f)", data$ebgm_ci_low, data$ebgm_ci_high),
      stringsAsFactors = FALSE
    )
  )
  
  if (method %in% names(column_definitions)) {
    return(column_definitions[[method]])
  } else {
    return(data.frame())  # Return empty data frame for unknown methods
  }
}

add_metadata_attributes <- function(table, term_col, method_config, 
                                    show_significant_only, clear_format) {
  
  attr(table, "term_type") <- term_col
  attr(table, "methods_used") <- method_config$used_methods
  attr(table, "significance_filter") <- show_significant_only
  attr(table, "format_type") <- ifelse(clear_format, "clear", "original")
  attr(table, "timestamp") <- Sys.time()
  
  return(table)
}
