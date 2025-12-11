#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]

# -------------------------------------------------------------------------
# Helper: safe read
# -------------------------------------------------------------------------
safe_read <- function(path) {
  if (length(path) == 0) return(NULL)
  df <- try(read_csv(path, show_col_types = FALSE), silent = TRUE)
  if (inherits(df, "try-error")) return(NULL)
  if (nrow(df) == 0) return(NULL)
  df
}

# -------------------------------------------------------------------------
# Locate files
# -------------------------------------------------------------------------
m1_file <- list.files(file.path(outdir, "processed_data"),
                      pattern = "dauerMod_NonOverlappingWorms.csv",
                      full.names = TRUE)

m2_file <- list.files(file.path(outdir, "processed_data"),
                      pattern = "nondauerMod_NonOverlappingWorms.csv",
                      full.names = TRUE)

# -------------------------------------------------------------------------
# Load files safely
# -------------------------------------------------------------------------
df1 <- safe_read(m1_file)
df2 <- safe_read(m2_file)

# If both empty, create placeholder
if (is.null(df1) && is.null(df2)) {
  message("WARNING: No dauer outputs found. Creating empty summary.")
  write_csv(tibble(Message="NO_WORMS_DETECTED"), 
            file.path(outdir, "dauer_summary.csv"))
  quit(save="no")
}

# Add model column
if (!is.null(df1)) df1 <- df1 %>% mutate(model="dauer")
if (!is.null(df2)) df2 <- df2 %>% mutate(model="nondauer")

# Merge available data
combined <- bind_rows(df1, df2)

# -------------------------------------------------------------------------
# Fix missing columns safely
# -------------------------------------------------------------------------
# If column missing, add NA
needed <- c("ImageNumber","Metadata_Date","Metadata_Well","Metadata_Plate")

for (nm in needed) {
  if (!nm %in% names(combined)) {
    combined[[nm]] <- NA
  }
}

# Clean Metadata_Date
combined <- combined %>%
  mutate(Metadata_Date = suppressWarnings(as.integer(Metadata_Date)))

# -------------------------------------------------------------------------
# Save outputs
# -------------------------------------------------------------------------
write_csv(combined, file.path(outdir, "dauer_all_raw.csv"))

summary_tbl <- combined %>%
  group_by(model) %>%
  summarise(count = n(), .groups="drop")

write_csv(summary_tbl, file.path(outdir, "dauer_summary.csv"))

message("✓ proc_CP_output_dauer.R completed successfully")
