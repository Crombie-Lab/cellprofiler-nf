#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
===========================================================
  VERSION CHECK
===========================================================
*/
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow >= 20.0 -- you are running $nextflow.version"
    System.exit(1)
}

/*
===========================================================
  PARAMETERS & DEFAULTS
===========================================================
*/
date = new Date().format('yyyyMMdd')

params.pipeline          = params.pipeline ?: null
params.project           = params.project  ?: "${workflow.projectDir}/input_data/20220420_testproject"
params.groups            = params.groups   ?: "plate,well"
params.mask              = params.mask     ?: "wellmask_98.png"
params.out               = params.out      ?: "${workflow.projectDir}/Analysis-${date}"

params.slurm_partition   = params.slurm_partition ?: "long"
params.array_size        = (params.array_size        ?: 6)  as int
params.array_concurrency = (params.array_concurrency ?: 3)  as int
params.cp_extra_args     = params.cp_extra_args ?: ""

params.data_dir       = "${workflow.projectDir}/input_data"
params.bin_dir        = "${workflow.projectDir}/bin"
params.raw_pipe_dir   = "${params.data_dir}/CP_pipelines"
params.well_mask_path = "${params.data_dir}/well_masks/${params.mask}"
params.worm_model_dir = "${params.data_dir}/worm_models"

/*
===========================================================
  PIPELINE MODE (TWO MODELS ONLY)
===========================================================
*/
def pipe
def worm_model1
def worm_model2
def model_name1
def model_name2
def config_script_name
def proc_script_name

switch (params.pipeline) {

    case 'dauer':
        pipe = "dauer-nf"
        worm_model1 = "dauerMod.xml"
        worm_model2 = "nondauerMod.xml"
        model_name1 = "dauerMod_NonOverlappingWorms"
        model_name2 = "nondauerMod_NonOverlappingWorms"
        config_script_name = "config_CP_input_dauer.R"
        proc_script_name   = "proc_CP_output_dauer.R"
        break

    case 'toxin':
        pipe = "toxin-nf-keyence"
        worm_model1 = "L4_N2_HB101_100w.xml"
        worm_model2 = "L2L3_N2_HB101_100w.xml"
        model_name1 = "L4_N2_HB101_100w_NonOverlappingWorms"
        model_name2 = "L2L3_N2_HB101_100w_NonOverlappingWorms"
        config_script_name = "config_CP_input_toxin.R"
        proc_script_name   = "proc_CP_output_toxin.R"
        break

    default:
        println "ERROR: specify --pipeline dauer OR --pipeline toxin"
        System.exit(1)
}

/*
===========================================================
  LOG INFO
===========================================================
*/
log.info """
C E L L P R O F I L E R  –  N E X T F L O W
===========================================

Pipeline         : ${params.pipeline}
Project dir      : ${params.project}
Groups           : ${params.groups}
Mask             : ${params.mask}
Out dir          : ${params.out}
SLURM partition  : ${params.slurm_partition}
Array size       : ${params.array_size}
Array concurrency: ${params.array_concurrency}
"""

/*
===========================================================
  CONFIG PROCESS (COMMON FOR BOTH PIPELINES)
===========================================================
*/
process config_CP_input {

    publishDir "${params.out}/pipeline", mode: 'copy', pattern: "pipeline.cppipe"
    publishDir "${params.out}/metadata", mode: 'copy', pattern: "metadata.csv"
    publishDir "${params.out}/groups",   mode: 'copy', pattern: "groups.tsv"

    input:
    tuple file(raw_pipe),
          val(project_dir),
          val(mask_path),
          val(group_cols),
          val(out_dir),
          val(worm_model_dir),
          val(model1),
          val(model2),
          val(config_mode),
          file(config_script)

    output:
    tuple path("groups.tsv"), val(out_dir), path("pipeline.cppipe"), emit: config_out
    path "metadata.csv", emit: metadata_file

    script:
    """
    set -Eeuo pipefail

    # Edit .cppipe
    awk '{gsub(/METADATA_DIR/,   "${out_dir}/metadata"); print}' ${raw_pipe} |
    awk '{gsub(/METADATA_CSV_FILE/, "metadata.csv"); print}' |
    awk '{gsub(/WORM_MODEL_DIR/, "${worm_model_dir}"); print}' |
    awk '{gsub(/MODEL1_XML_FILE/, "${model1}"); print}' |
    awk '{gsub(/MODEL2_XML_FILE/, "${model2}"); print}' \
        > pipeline.cppipe

    if [ "${config_mode}" = "toxin" ]; then
        tmp=\$(mktemp)
        awk '{gsub(/WELL_MASK/, "${mask_path}"); print}' pipeline.cppipe > \$tmp
        mv \$tmp pipeline.cppipe
    fi

    # Run R script
    Rscript --vanilla ${config_script} \
        "${project_dir}" \
        "${mask_path}" \
        "${group_cols}" \
        "${out_dir}/pipeline/pipeline.cppipe" \
        "${out_dir}"
    """
}

/*
===========================================================
  SLURM ARRAY EXECUTION (ONE JOB, MANY TASKS)
===========================================================
*/
process CP_Array {

    executor 'slurm'
    queue params.slurm_partition

    cpus 2
    memory '16 GB'
    time '4h'

    clusterOptions "--array=1-${params.array_size}%${params.array_concurrency}"

    input:
    tuple path(groups_file), val(out_dir), path(pipeline_file)

    output:
    tuple val("done"), val(out_dir), emit: cp_out

    script:
    """
    set -Eeuo pipefail

    i=\${SLURM_ARRAY_TASK_ID:-1}

    maxrows=\$(awk 'END{print NR-1}' "${groups_file}")

    if (( i < 1 || i > maxrows )); then
        exit 0
    fi

    line=\$(awk -v idx="\$i" 'NR==1{next} NR==idx+1{print; exit}' "${groups_file}")

    group=\$(echo "\$line"  | awk -F '\\t' '{print \$1}')
    output=\$(echo "\$line" | awk -F '\\t' '{print \$3}')

    mkdir -p "\$output"

    cellprofiler -c -r -p "${pipeline_file}" -g "\$group" -o "\$output" ${params.cp_extra_args}
    """
}

/*
===========================================================
  MERGE OUTPUTS — DAUER
===========================================================
*/
process proc_CP_output_dauer {

    input:
    tuple val(done), val(out_dir),
          val(m1), val(m2),
          file(proc_script)

    script:
    """
    set -Eeuo pipefail

    mkdir -p "${out_dir}/processed_data"
    mkdir -p "${out_dir}/processed_images"

    find "${out_dir}/CP_output" -type f -name "${m1}.csv" -print0 |
        xargs -0 awk 'FNR>1 || NR==1 {print}' \
        > "${out_dir}/processed_data/${m1}.csv" || true

    find "${out_dir}/CP_output" -type f -name "${m2}.csv" -print0 |
        xargs -0 awk 'FNR>1 || NR==1 {print}' \
        > "${out_dir}/processed_data/${m2}.csv" || true

    find "${out_dir}/CP_output" -type f -name "*.png" -print0 |
        xargs -0 -I{} mv "{}" "${out_dir}/processed_images" || true

    Rscript --vanilla "${proc_script}" "${out_dir}"
    """
}

/*
===========================================================
  MERGE OUTPUTS — TOXIN
===========================================================
*/
process proc_CP_output_toxin {

    input:
    tuple val(done), val(out_dir),
          val(m1), val(m2),
          file(proc_script)

    script:
    """
    set -Eeuo pipefail

    mkdir -p "${out_dir}/processed_data"
    mkdir -p "${out_dir}/processed_images"

    find "${out_dir}/CP_output" -type f -name "${m1}.csv" -print0 |
        xargs -0 awk 'FNR>1 || NR==1 {print}' \
        > "${out_dir}/processed_data/${m1}.csv" || true

    find "${out_dir}/CP_output" -type f -name "${m2}.csv" -print0 |
        xargs -0 awk 'FNR>1 || NR==1 {print}' \
        > "${out_dir}/processed_data/${m2}.csv" || true

    find "${out_dir}/CP_output" -type f -name "WormObjects.csv" -print0 |
        xargs -0 awk 'FNR>1 || NR==1 {print}' \
        > "${out_dir}/processed_data/WormObjects.csv" || true

    find "${out_dir}/CP_output" -type f -name "*.png" -print0 |
        xargs -0 -I{} mv "{}" "${out_dir}/processed_images" || true

    Rscript --vanilla "${proc_script}" "${out_dir}"
    """
}

/*
===========================================================
  WORKFLOW
===========================================================
*/
workflow {

    config_input_ch = Channel
        .fromPath("${params.raw_pipe_dir}/${pipe}.cppipe")
        .combine(Channel.of(params.project))
        .combine(Channel.of(params.well_mask_path))
        .combine(Channel.of(params.groups))
        .combine(Channel.of(params.out))
        .combine(Channel.of(params.worm_model_dir))
        .combine(Channel.of(worm_model1))
        .combine(Channel.of(worm_model2))
        .combine(Channel.of(params.pipeline))
        .combine(Channel.fromPath("${params.bin_dir}/${config_script_name}"))

    configured = config_CP_input(config_input_ch)

    cp_results = CP_Array(configured.config_out)

    if (params.pipeline == 'dauer') {

        proc_input = cp_results.cp_out
            .combine(Channel.of(model_name1))
            .combine(Channel.of(model_name2))
            .combine(Channel.fromPath("${params.bin_dir}/${proc_script_name}"))

        proc_CP_output_dauer(proc_input)
    }
    else {

        proc_input = cp_results.cp_out
            .combine(Channel.of(model_name1))
            .combine(Channel.of(model_name2))
            .combine(Channel.fromPath("${params.bin_dir}/${proc_script_name}"))

        proc_CP_output_toxin(proc_input)
    }
}

/*
===========================================================
  SUMMARY
===========================================================
*/
workflow.onComplete {
    println """
    Pipeline Summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration     : ${workflow.duration}
    Success      : ${workflow.success}
    WorkDir      : ${workflow.workDir}
    """
}
