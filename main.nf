#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
===========================================================
  PARAMETERS — HARD DEFAULTS (NO WARNINGS)
===========================================================
*/

date = new Date().format('yyyyMMdd')

params.pipeline          = params.pipeline ?: null
params.project           = params.project  ?: null
params.groups            = params.groups   ?: "plate,well"
params.mask              = params.mask     ?: "wellmask_98.png"

params.slurm_partition   = params.slurm_partition   ?: "long"
params.cp_extra_args     = params.cp_extra_args     ?: ""

params.data_dir  = "${workflow.projectDir}/input_data"
params.bin_dir   = "${workflow.projectDir}/bin"
params.pipe_dir  = "${params.data_dir}/CP_pipelines"
params.model_dir = "${params.data_dir}/worm_models"
params.mask_path = "${params.data_dir}/well_masks/${params.mask}"

if( !params.pipeline || !params.project ) {
    error "You MUST specify --pipeline and --project"
}

params.out = "${params.project}/Analysis-${date}"

/*
===========================================================
  PIPELINE SELECTION
===========================================================
*/

def pipe
def config_script
def proc_script

if (params.pipeline == "toxin") {

    pipe          = "toxin-nf-keyence.cppipe"
    config_script = "config_CP_input_toxin.R"
    proc_script   = "proc_CP_output_toxin.R"

}
else if (params.pipeline == "dauer") {

    pipe          = "dauer-nf.cppipe"
    config_script = "config_CP_input_dauer.R"
    proc_script   = "proc_CP_output_dauer.R"

}
else {
    error "Invalid --pipeline (use toxin or dauer)"
}

/*
===========================================================
  BANNER (PRINT AFTER PARAMS EXIST)
===========================================================
*/

println """
===========================================================
 CELLPROFILER - NEXTFLOW
===========================================================
Pipeline        : ${params.pipeline}
Project dir     : ${params.project}
Groups          : ${params.groups}
Mask            : ${params.mask}
Output dir      : ${params.out}
SLURM partition : ${params.slurm_partition}
Array size      : ${params.array_size}
Array conc.     : ${params.array_concurrency}
===========================================================
"""

/*
===========================================================
  CONFIGURATION
===========================================================
*/

process CONFIG_CP {

    publishDir "${params.out}/pipeline", mode: 'copy'
    publishDir "${params.out}/metadata", mode: 'copy'
    publishDir "${params.out}/groups",   mode: 'copy'

    input:
    file raw_pipe
    file config_r
    val  project
    val  groups
    val  out_dir

    output:
    tuple path("groups.tsv"), val(out_dir), path("pipeline.cppipe")

    script:
    """
    awk '{gsub(/METADATA_DIR/,"${out_dir}/metadata"); print}' ${raw_pipe} |
    awk '{gsub(/METADATA_CSV_FILE/,"metadata.csv"); print}' |
    awk '{gsub(/WORM_MODEL_DIR/,"${params.model_dir}"); print}' |
    awk '{gsub(/WELL_MASK/,"${params.mask}"); print}' > pipeline.cppipe

    Rscript --vanilla ${config_r} \
        "${project}" \
        "${params.mask_path}" \
        "${groups}" \
        pipeline.cppipe \
        "${out_dir}"
    """
}

/*
===========================================================
  SLURM ARRAY EXECUTION
===========================================================
*/

process CP_ARRAY {

    executor 'slurm'
    queue params.slurm_partition
    clusterOptions "--array=1-${params.array_size}%${params.array_concurrency}"

    input:
    tuple path(groups_tsv), val(out_dir), path(pipeline)

    output:
    val(out_dir)

    script:
    """
    i=\${SLURM_ARRAY_TASK_ID:-1}

    line=\$(awk -v i="\$i" 'NR==i+1{print}' ${groups_tsv})
    [ -z "\$line" ] && exit 0

    group=\$(echo "\$line" | cut -f1)
    output=\$(echo "\$line" | cut -f3)

    mkdir -p "\$output"

    cellprofiler -c -r \
        -p pipeline.cppipe \
        -g "\$group" \
        -o "\$output"
    """
}

/*
===========================================================
  POST PROCESSING
===========================================================
*/

process POST {

    input:
    val out_dir

    script:
    """
    Rscript --vanilla ${params.bin_dir}/${proc_script} "${out_dir}"
    """
}

/*
===========================================================
  WORKFLOW
===========================================================
*/

workflow {

    cfg = CONFIG_CP(
        file("${params.pipe_dir}/${pipe}"),
        file("${params.bin_dir}/${config_script}"),
        params.project,
        params.groups,
        params.out
    )

    run = CP_ARRAY(cfg)

    POST(run)
}

workflow.onComplete {
    println """
===========================================================
 Pipeline Summary
===========================================================
Completed at : ${workflow.complete}
Duration     : ${workflow.duration}
Success      : ${workflow.success}
WorkDir      : ${workflow.workDir}
===========================================================
"""
}
