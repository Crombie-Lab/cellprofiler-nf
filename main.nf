nextflow.enable.dsl=2

// ---------- params (safe defaults)
params.pipe              = params.pipe               ?: "$HOME/cellprofiler-nf/input_data/20220420_testproject/Analysis-20251110/pipeline/pipeline.cppipe"
params.outdir            = params.outdir             ?: "$HOME/cellprofiler-nf/input_data/20220420_testproject/out"
params.array_size        = (params.array_size        ?: 100) as int
params.array_concurrency = (params.array_concurrency ?: 10)  as int
params.slurm_partition   = params.slurm_partition    ?: 'med'
params.cp_extra_args     = params.cp_extra_args      ?: ''   // silence the WARN

// ---------- process: one NF task -> one SLURM job array
process CP_Array {
  tag "slurm-array"
  executor 'slurm'
  queue params.slurm_partition
  cpus 1
  memory '2 GB'
  time '2h'

  // one sbatch with N children
  clusterOptions { "--array=1-${params.array_size}%${params.array_concurrency}" }

  // make CellProfiler available on compute nodes
  beforeScript """
    set -Eeuo pipefail
    if [ -f "\$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
      . "\$HOME/miniforge3/etc/profile.d/conda.sh" || true
      conda activate cp-426 || true
      export LD_LIBRARY_PATH="\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH"
    fi
  """

  // no inputs; we read params.* directly
  script:
  """
  set -Eeuo pipefail

  i="\${SLURM_ARRAY_TASK_ID:-1}"
  OUTDIR="${params.outdir}/\${i}"
  PIPE="${params.pipe}"

  [[ -f "\$PIPE" ]] || { echo "[FATAL] Missing pipeline: \$PIPE"; exit 2; }
  mkdir -p "\$OUTDIR"

  echo "[INFO] Task=\$i"
  echo "[INFO] Pipeline: \$PIPE"
  echo "[INFO] Output  : \$OUTDIR"

  command -v cellprofiler >/dev/null 2>&1 || { echo "[FATAL] cellprofiler not found on PATH"; exit 127; }

  cellprofiler -c -r -p "\$PIPE" -o "\$OUTDIR" ${params.cp_extra_args}

  echo "[DONE] Task \$i -> \$OUTDIR"
  """
}
