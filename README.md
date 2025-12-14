Cellprofiler-nf

A Nextflow pipeline to run CellProfiler pipelines on raw images with SLURM array execution and standardized post-processing.

Pipeline overview

cellprofiler-nf is a Nextflow + CellProfiler pipeline designed to process large microscopy image datasets using SLURM array jobs.
Each metadata group (e.g., plate/well) is processed as a separate array task, allowing efficient parallel execution on HPC systems such as ai-panther.
The pipeline supports two CellProfiler workflows:
-	dauer
-	toxin
Both pipelines:
-	Read raw images
-	Generate metadata and grouping
-	Run CellProfiler headless
-	Merge outputs into clean, analysis-ready results
---
Key features
-	SLURM array job execution (one task per metadata group)
-	Supports dauer and toxin CellProfiler pipelines
-	Clean, reproducible directory structure
-	Automatic metadata + grouping generation
-	Post-processing of CSV outputs and images
-	Designed for HPC use, not laptop hacks

---
ai-panther usage
# clone the repository
git clone https://github.com/AndersenLab/cellprofiler-nf.git
cd cellprofiler-nf

# load required modules
module load singularity

# activate your Nextflow conda environment
source activate nextflow-23.10
Example runs
Run the toxin pipeline
nextflow run main.nf \
  --pipeline toxin \
  --project /path/to/20220501_toxinDebug
Run the dauer pipeline
nextflow run main.nf \
  --pipeline dauer \
  --project /path/to/20220501_dauerDebug
Override SLURM array settings at runtime
nextflow run main.nf \
  --pipeline dauer \
  --project /path/to/project \
  --array_size 6 \
  --array_concurrency 3
________________________________________
cellprofiler-nf help
C E L L P R O F I L E R  -  N E X T F L O W
===============================================

Usage:
nextflow run main.nf --pipeline <dauer|toxin> --project <project directory>

Mandatory arguments:
--pipeline            CellProfiler pipeline to use: dauer or toxin
--project             Path to the project directory

Optional arguments:
--groups              Metadata grouping (default: plate,well)
--mask                Well mask filename (default: wellmask_98.png)
--out                 Output directory (default: project/Analysis-YYYYMMDD)
--array_size          SLURM array size (default set in main.nf)
--array_concurrency   Max concurrent array tasks
--slurm_partition     SLURM partition (default: long)
________________________________________
Input directory structure
Each project directory must contain a raw_images/ folder.
Dauer input example
20220501_dauerDebug/
├── raw_images
│   ├── 20220501_dauerDebug-p002-m2X_A01_w1.TIF
│   ├── 20220501_dauerDebug-p002-m2X_A01_w2.TIF
│   └── ...
Toxin input example
20220501_toxinDebug/
├── raw_images
│   ├── 20220501_toxinDebug-p001-m2X_A01.TIF
│   ├── 20220501_toxinDebug-p010-m2X_A01.TIF
│   └── ...
Naming conventions
•	Dauer
Date-Experiment-Plate-Magnification_Well_Wavelength.TIF
•	Toxin
Date-Experiment-Plate-Magnification_Well.TIF
---
Output directory structure
By default, outputs are written to:
<project>/Analysis-YYYYMMDD/
Common output layout
Analysis-YYYYMMDD/
├── pipeline/
│   └── pipeline.cppipe
├── metadata/
│   └── metadata.csv
├── groups/
│   └── groups.tsv
├── CP_output/
│   └── <per-group CellProfiler outputs>
├── processed_data/
├── processed_images/
---_
Dauer output example
processed_data/
└── 20220501_dauerDebug_Analysis-YYYYMMDD.RData

processed_images/
├── *_overlay.png
├── *_dauerMod_straightened_RFP.png
├── *_nondauerMod_straightened_RFP.png
└── ...
________________________________________
Toxin output example
processed_data/
└── 20220501_toxinDebug_Analysis-YYYYMMDD.RData

processed_images/
├── *_overlay.png
└── ...
________________________________________
How SLURM array execution works
1.	config_CP_input:
o	Builds pipeline.cppipe
o	Generates metadata.csv
o	Generates groups.tsv
2.	CP_Array:
o	Submits one SLURM job
o	Each array task processes one row of groups.tsv
o	Runs CellProfiler headless per group
3.	proc_CP_output_*:
o	Merges CSV outputs
o	Collects images
o	Generates final RData files
This design scales to hundreds of wells without code changes.
________________________________________
Dependencies
Option 1: HPC environment (recommended)
-	Nextflow ≥ 20
-	SLURM
-	Singularity
-	R (with required packages)
-	CellProfiler (inside Singularity image)
On ai-panther, these are already available.
________________________________________
Option 2: Docker (local testing)
If using Docker, you do not need to install dependencies manually.
Example (conceptual):
docker run -it \
  -v $PWD:/workspace \
  cellprofiler/cellprofiler:latest \
  nextflow run main.nf --pipeline dauer --project /workspace/example
