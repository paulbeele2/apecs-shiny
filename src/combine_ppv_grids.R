#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

base_dir <- "results/rshiny_slider"
output_als <- "./varying_param_als_grid.csv"
output_alsftd <- "./varying_param_alsftd_grid.csv"

if (!dir.exists(base_dir)) {
  stop("Base directory does not exist: ", base_dir)
}

grid_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

if (length(grid_dirs) == 0) {
  stop("No parameter directories found in: ", base_dir)
}

extract_params <- function(dir_path) {
  dir_name <- basename(dir_path)

  m <- str_match(
    dir_name,
    "^ALSFTD_common_([0-9.]+)_ALSFTD_rare_([0-9.]+)_h2als_([0-9.]+)$"
  )

  if (any(is.na(m))) {
    stop("Directory name does not match expected pattern: ", dir_name)
  }

  tibble(
    ALSFTD_common = as.numeric(m[2]),
    ALSFTD_rare   = as.numeric(m[3]),
    h2als         = as.numeric(m[4]),
    source_dir    = dir_path
  )
}

read_grid_with_params <- function(dir_path, filename) {
  file_path <- file.path(dir_path, filename)

  if (!file.exists(file_path)) {
    warning("Missing file, skipping: ", file_path)
    return(NULL)
  }

  params <- extract_params(dir_path)

  df <- read_csv(file_path, show_col_types = FALSE)

  bind_cols(params %>% select(-source_dir), df)
}

als_list <- map(grid_dirs, read_grid_with_params, filename = "output_ppv_als_grid.csv")
alsftd_list <- map(grid_dirs, read_grid_with_params, filename = "output_ppv_alsftd_grid.csv")

als_combined <- bind_rows(als_list) %>%
  arrange(
    ALSFTD_common, ALSFTD_rare, h2als,
    relatives_1st_als, relatives_2nd_als, relatives_3rd_als
  )

alsftd_combined <- bind_rows(alsftd_list) %>%
  arrange(
    ALSFTD_common, ALSFTD_rare, h2als,
    relatives_1st_als, relatives_2nd_als, relatives_3rd_als,
    relatives_1st_ftd_unique, relatives_2nd_ftd_unique, relatives_3rd_ftd_unique
  )

write_csv(als_combined, output_als)
write_csv(alsftd_combined, output_alsftd)

message("Wrote ALS grid to: ", output_als, " (rows: ", nrow(als_combined), ")")
message("Wrote ALS+FTD grid to: ", output_alsftd, " (rows: ", nrow(alsftd_combined), ")")