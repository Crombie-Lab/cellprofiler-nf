#!/usr/bin/env bash
# Disable 'set -u' until after /etc/profile to avoid LC_* errors
set -e
source /etc/profile
module load apptainer/1.3.4-gcc-14.2.0-spxhran
set -u  # Re-enable after module load

PROJ="$HOME/cellprofiler-nf/input_data/20220420_testproject"
PIPE="$PROJ/Analysis-20251110/pipeline/pipeline.cppipe"
ID="${SLURM_ARRAY_TASK_ID:?missing array index}"
ANALYSIS_DIR="$PROJ/Analysis-array-$ID"
EXP_CSV="$ANALYSIS_DIR/Experiment.csv"
OUT="$ANALYSIS_DIR/_cp_out"
REP="$PROJ/reports"
LOG="$REP/trace-$ID.txt"
SIF="/home1/hmohammedshe2021/cp-4.2.6.sif"

mkdir -p "$OUT" "$REP"

if [ ! -s "$EXP_CSV" ]; then
  echo "FATAL: $EXP_CSV missing or empty" | tee "$LOG"
  exit 2
fi

echo "[$(date)] START id=$ID" | tee "$LOG"

apptainer exec --no-home -B "$HOME" --cleanenv "$SIF" \
  cellprofiler -c -r -p "$PIPE" -o "$OUT" --data-file "$EXP_CSV" \
  2>&1 | tee -a "$LOG"

printf "<html><body><h1>Array %s OK</h1></body></html>\n" "$ID" > "$REP/array-$ID.html"
echo "[$(date)] DONE id=$ID" | tee -a "$LOG"
