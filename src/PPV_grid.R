#!/usr/bin/env Rscript

## PPV_grid.R
## Usage:
## Rscript scripts/PPV_grid.R ../results/PPV_grid/combined_simulations.csv \
##                    output_ppv_als_grid.csv \
##                    output_ppv_alsftd_grid.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript scripts/PPV_grid.R <input_combined_simulations.csv> <output_ppv_als_grid.csv> <output_ppv_alsftd_grid.csv>")
}

input_file         <- args[1]
output_als_file    <- args[2]
output_alsftd_file <- args[3]

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

message("Reading combined simulations from: ", input_file)
results_df <- read_csv(input_file, show_col_types = FALSE)

## Wilson CI exactly as in prediction_model.R
wilson_ci <- function(x, n, conf_level = 0.95) {
  if (n == 0 || x < 0 || x > n) return(c(NA, NA, NA))
  
  z  <- qnorm(1 - (1 - conf_level) / 2)
  z2 <- z^2
  p_hat <- x / n

  center <- (x + z2 / 2) / (n + z2)
  margin <- z * sqrt((x * (n - x) + z2 / 4) / (n + z2)^2) / sqrt(n)
  
  ci_lower <- center - margin
  ci_upper <- center + margin
  
  ci_lower <- pmax(0, ci_lower)
  ci_upper <- pmin(1, ci_upper)
  
  c(estimate = p_hat, ci_lower = ci_lower, ci_upper = ci_upper)
}

calculate_ppv_for_combinations <- function(results_df, vars, max_count = 4, conf_level = 0.95) {
  missing_vars <- setdiff(c(vars, "mendel_ALS_Y"), names(results_df))
  if (length(missing_vars) > 0) {
    stop("Missing required columns in input: ",
         paste(missing_vars, collapse = ", "))
  }
  
  df <- results_df[, c(vars, "mendel_ALS_Y")]
  total_n <- nrow(df)
  
  # Cap counts
  for (v in vars) {
    df[[v]] <- pmin(df[[v]], max_count)
  }
  
  combo_counts <- df %>%
    group_by(across(all_of(vars))) %>%
    summarise(
      n = n(),
      n_mendelian = sum(mendel_ALS_Y == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      n_non_mendelian = n - n_mendelian
    )
  
  ppv_mat <- t(vapply(
    seq_len(nrow(combo_counts)),
    function(i) wilson_ci(combo_counts$n_mendelian[i],
                          combo_counts$n[i],
                          conf_level = conf_level),
    FUN.VALUE = c(estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
  ))
  
  # Human-readable labels, matching your earlier pattern
  criterion_label <- with(combo_counts, {
    if ("relatives_1st_ftd_unique" %in% vars) {
      paste0(
        "1st_ALS=", relatives_1st_als,
        "; 2nd_ALS=", relatives_2nd_als,
        "; 3rd_ALS=", relatives_3rd_als,
        "; 1st_FTD=", relatives_1st_ftd_unique,
        "; 2nd_FTD=", relatives_2nd_ftd_unique,
        "; 3rd_FTD=", relatives_3rd_ftd_unique
      )
    } else {
      paste0(
        "1st_ALS=", relatives_1st_als,
        "; 2nd_ALS=", relatives_2nd_als,
        "; 3rd_ALS=", relatives_3rd_als
      )
    }
  })
  
  combo_counts %>%
    mutate(
      PPV         = round(ppv_mat[, 1], 3),
      PPV_CI_low  = round(ppv_mat[, 2], 3),
      PPV_CI_high = round(ppv_mat[, 3], 3),
      prevalence  = sprintf("%d out of %d simulated index patients", n, total_n),
      criterion_label = criterion_label
    ) %>%
    arrange(
      relatives_1st_als,
      relatives_2nd_als,
      relatives_3rd_als
    )
}

## ALS-only grid
als_vars <- c(
  "relatives_1st_als",
  "relatives_2nd_als",
  "relatives_3rd_als"
)

message("Computing PPV grid (ALS only)...")
ppv_als_grid <- calculate_ppv_for_combinations(
  results_df,
  vars = als_vars,
  max_count = 5,
  conf_level = 0.95
)

message("Writing ALS PPV grid to: ", output_als_file)
write_csv(ppv_als_grid, output_als_file)
message("ALS grid rows: ", nrow(ppv_als_grid))

## ALS + FTD grid
alsftd_vars <- c(
  "relatives_1st_als",
  "relatives_2nd_als",
  "relatives_3rd_als",
  "relatives_1st_ftd_unique",
  "relatives_2nd_ftd_unique",
  "relatives_3rd_ftd_unique"
)

message("Computing PPV grid (ALS + FTD)...")
ppv_alsftd_grid <- calculate_ppv_for_combinations(
  results_df,
  vars = alsftd_vars,
  max_count = 5,
  conf_level = 0.95
)

message("Writing ALS+FTD PPV grid to: ", output_alsftd_file)
write_csv(ppv_alsftd_grid, output_alsftd_file)
message("ALS+FTD grid rows: ", nrow(ppv_alsftd_grid))

message("Done.")