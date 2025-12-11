#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(glue)
  library(purrr)
  library(data.table)
})

# Args:
# 1 - full path to project directory (with raw_images/)
# 2 - full path to well mask (e.g. .../well_masks/wellmask_98.png)
# 3 - group cols string, e.g. "plate,well"
# 4 - edited pipeline path (pipeline.cppipe)  [used only for groups.tsv]
# 5 - out path (e.g. Analysis-YYYYMMDD)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: config_CP_input_dauer.R <project_dir> <mask_path> <groups> <pipeline_path> <out_dir>")
}

projDir       <- args[1]
mask_path     <- args[2]
group_arg     <- args[3]
pipeline_path <- args[4]
out_dir       <- args[5]

message("Project dir  : ", projDir)
message("Mask path    : ", mask_path)
message("Groups       : ", group_arg)
message("Pipeline path: ", pipeline_path)
message("Out dir      : ", out_dir)

raw_imagesDir <- file.path(projDir, "raw_images")

if (!dir.exists(raw_imagesDir)) {
  stop("raw_images directory does not exist: ", raw_imagesDir)
}

files <- list.files(path = raw_imagesDir, full.names = FALSE)
if (length(files) == 0) {
  stop("No image files found in raw_images: ", raw_imagesDir)
}

meta1 <- tibble(
  file      = list.files(path = raw_imagesDir, full.names = FALSE),
  file_path = list.files(path = raw_imagesDir, full.names = TRUE)
) %>%
  mutate(copy = file) %>%
  # e.g. 20220420-testproject-p002-m2X_A01_w1.TIF
  separate(col = copy, into = c("date", "exp", "plate", "mag"), sep = "-") %>%
  separate(col = mag,   into = c("mag", "well", "wave"), sep = "_") %>%
  separate(col = wave,  into = c("wave", "TIF"), sep = "[.]") %>%
  select(-TIF) %>%
  mutate(
    row = str_extract(well, "[A-Z]"),
    col = str_extract(well, "[0-9][0-9]"),
    Image_PathName_wellmask_98.png = str_replace(mask_path, "([^/]+$)", ""),
    Image_FileName_wellmask_98.png = str_extract(mask_path, "([^/]+$)")
  )

# how many wavelengths (w1, w2, ...)
n_wave <- length(unique(meta1$wave))
message("Detected ", n_wave, " wavelength(s): ", paste(unique(meta1$wave), collapse = ", "))

# groups: e.g. "plate,well"
groups <- str_split(group_arg, pattern = ",")[[1]]
meta1$group <- apply(meta1[, groups, drop = FALSE], 1, paste, collapse = "_")

meta2 <- meta1 %>%
  pivot_wider(names_from = wave, values_from = c(file, file_path)) %>%
  # assume w1 = BF, w2 = RFP (matches your dataset)
  rename(
    Image_FileName_RawBF = file_w1,
    Image_PathName_RawBF = file_path_w1,
    Image_FileName_RawRFP = file_w2,
    Image_PathName_RawRFP = file_path_w2
  ) %>%
  mutate(
    Image_PathName_RawRFP = str_replace(Image_PathName_RawRFP, "([^/]+$)", ""),
    Image_PathName_RawBF  = str_replace(Image_PathName_RawBF,  "([^/]+$)", "")
  ) %>%
  select(
    Metadata_Experiment = exp,
    Metadata_Date       = date,
    Metadata_Plate      = plate,
    Metadata_Well       = well,
    Metadata_Group      = group,
    Metadata_Magnification = mag,
    Image_FileName_RawBF,
    Image_PathName_RawBF,
    Image_FileName_RawRFP,
    Image_PathName_RawRFP,
    Image_FileName_wellmask_98.png,
    Image_PathName_wellmask_98.png
  )

# Write metadata.csv in current working dir (Nextflow publishDir will copy it)
write.table(
  meta2,
  file      = "metadata.csv",
  quote     = FALSE,
  sep       = ",",
  row.names = FALSE
)

# groups.tsv: one row per group, with CP_output path under out_dir
gs <- meta2 %>%
  distinct(Metadata_Group, .keep_all = TRUE) %>%
  mutate(
    group   = paste0("Metadata_Group=", Metadata_Group),
    pipeline = pipeline_path,
    output  = paste0(out_dir, "/CP_output/", Metadata_Group)
  ) %>%
  select(group:output)

write.table(
  gs,
  file      = "groups.tsv",
  quote     = FALSE,
  sep       = "\t",
  row.names = FALSE
)

# Make CP_output dirs
for (i in unique(gs$output)) {
  dir.create(i, recursive = TRUE, showWarnings = FALSE)
}

cat("✓ Metadata + groups.tsv created successfully\n")
