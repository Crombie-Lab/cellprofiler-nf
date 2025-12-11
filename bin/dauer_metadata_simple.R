#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(glue)

args <- commandArgs(trailingOnly = TRUE)

projDir  <- args[1]
maskPath <- args[2]
groupDef <- args[3]
pipeline <- args[4]
outDir   <- args[5]

# ------------------------------------------------------------------
# 1. Detect raw image folder (works even if names differ)
# ------------------------------------------------------------------
rawDir1 <- file.path(projDir, "raw_images")
rawDir2 <- file.path(projDir, "images")
rawDir3 <- projDir

raw_imagesDir <- if (dir.exists(rawDir1)) rawDir1 else
                  if (dir.exists(rawDir2)) rawDir2 else rawDir3

# ------------------------------------------------------------------
# 2. Load all image filenames
# ------------------------------------------------------------------
files <- list.files(raw_imagesDir, full.names = TRUE)

if (length(files) == 0)
    stop("No image files found. Check folder structure.", call. = FALSE)

df <- tibble(
    file = basename(files),
    file_path = dirname(files)
)

# ------------------------------------------------------------------
# 3. Extract WELL and PLATE from filenames (fallback safe)
# ------------------------------------------------------------------
df <- df %>%
  mutate(
    plate = "P1",
    well  = str_extract(file, "[A-H][0-9]{2}"),
    well  = ifelse(is.na(well), "A01", well),
    group = well
  )

# ------------------------------------------------------------------
# 4. Add mask name columns required by CellProfiler
# ------------------------------------------------------------------
mask_name <- basename(maskPath)
mask_dir  <- dirname(maskPath)

df$Image_FileName_wellmask <- mask_name
df$Image_PathName_wellmask <- mask_dir

# ------------------------------------------------------------------
# 5. Rename for CellProfiler expected metadata
# ------------------------------------------------------------------
meta <- df %>%
  mutate(
    Metadata_Experiment  = "DAUER",
    Metadata_Date        = "2024",
    Metadata_Plate       = plate,
    Metadata_Well        = well,
    Metadata_Group       = group,
    Metadata_Magnification = "10X",
    Image_FileName_RawBF = file,
    Image_PathName_RawBF = file_path
  ) %>%
  select(
    Metadata_Experiment,
    Metadata_Date,
    Metadata_Plate,
    Metadata_Well,
    Metadata_Group,
    Metadata_Magnification,
    Image_FileName_RawBF,
    Image_PathName_RawBF,
    Image_FileName_wellmask,
    Image_PathName_wellmask
  )

write_csv(meta, file.path(projDir, "metadata.csv"))

# ------------------------------------------------------------------
# 6. Create groups.tsv
# ------------------------------------------------------------------
gs <- meta %>%
  distinct(Metadata_Group, .keep_all = TRUE) %>%
  mutate(
    group    = paste0("Metadata_Group=", Metadata_Group),
    pipeline = pipeline,
    output   = file.path(outDir, "CP_output", Metadata_Group)
  ) %>%
  select(group, pipeline, output)

write.table(gs, file = file.path(projDir, "groups.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------
# 7. Create output directories
# ------------------------------------------------------------------
for (d in gs$output) dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("✓ Metadata + groups.tsv created successfully")
