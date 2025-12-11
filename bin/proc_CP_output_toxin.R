#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]

cat("[INFO] R: Starting toxin aggregation\n")

csv_dir <- file.path(out_dir, "processed_data")
if (!dir.exists(csv_dir)) {
  stop("processed_data directory does not exist: ", csv_dir)
}

files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  cat("[WARN] No CSV files found in processed_data\n")
  quit(save="no")
}

cat("[INFO] R: Found", length(files), "CSV files\n")

merge_all <- NULL

for (f in files) {
  cat("[INFO] Reading:", f, "\n")
  df <- tryCatch(
    read.csv(f, stringsAsFactors = FALSE),
    error = function(e) {
      cat("[WARN] Failed to read", f, "\n")
      return(NULL)
    }
  )
  
  if (!is.null(df) && nrow(df) > 0) {
    df$SourceModel <- basename(f)
    merge_all <- rbind(merge_all, df)
  }
}

if (is.null(merge_all)) {
  cat("[WARN] No valid CSV data to merge\n")
  quit(save="no")
}

outfile <- file.path(out_dir, "ALL_MODELS_MERGED.csv")
write.csv(merge_all, outfile, row.names = FALSE)

cat("[INFO] R: Finished. Output written to:\n")
cat(outfile, "\n")
