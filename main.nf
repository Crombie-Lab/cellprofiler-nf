#!/usr/bin/env nextflow

// Use DSL2
nextflow.enable.dsl=2

// QUEST nextflow version message
if( !nextflow.version.matches('>20.0') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

/*
~ ~ ~ > * PARAMETERS SETUP
*/

// Variables
date = new Date().format('yyyyMMdd')
// model_name = "NonOverlappingWorms" // JUST FOR NOW

// Setup pipeline parameter
params.pipeline = null
if("${params.pipeline}" == "dauer") {
    pipe = "dauer-nf"
    worm_model1 = "dauerMod.xml"
    worm_model2 = "nondauerMod.xml"
    model_name1 = "dauerMod_NonOverlappingWorms"
    model_name2 = "nondauerMod_NonOverlappingWorms"
} else if( "${params.pipeline}" == "toxin" ) {
    pipe = "toxin-nf-keyence" // To process keyence images with edgecorr/edgestrong circa 20241206 use "toxin-nf-keyence_edgestrong"
    worm_model1 = "L4_N2_HB101_100w.xml"
    worm_model2 = "L2L3_N2_HB101_100w.xml"
    worm_model3 = "L1_N2_HB101_100w.xml"
    worm_model4 = "MDHD.xml"
    model_name1 = "L4_N2_HB101_100w_NonOverlappingWorms"
    model_name2 = "L2L3_N2_HB101_100w_NonOverlappingWorms"
    model_name3 = "L1_N2_HB101_100w_NonOverlappingWorms"
    model_name4 = "MDHD_NonOverlappingWorms"
} else if(!params.pipeline) {
    println """
            Error: pipeline parameter not specified. Please enter --pipeline dauer or --pipeline toxin in command.
            """
            System.exit(1)
} else if("${params.pipeline}" != "toxin" || "${params.pipeline}" != "dauer" ) {
    println """
            Error: pipeline (${params.pipeline}) does not match expected value. Please enter either dauer or toxin.
            """
            System.exit(1)
}

// Configure other parameters
params.help = null
params.debug = null
params.project = null
params.groups = "plate,well"
params.mask = "20250209_well_mask.png"
params.data_dir = "${workflow.projectDir}/input_data" // this is different for gcp
params.bin_dir = "${workflow.projectDir}/bin" // this is different for gcp
params.well_mask = "${params.data_dir}/well_masks/${params.mask}"
params.out = "${params.project}/Analysis-${date}"
params.raw_pipe_dir = "${params.data_dir}/CP_pipelines"
params.raw_pipe = "${params.raw_pipe_dir}/${pipe}.cppipe"
params.edited_pipe = "${params.out}/pipeline/pipeline.cppipe"
params.metadata_dir = "${params.out}/metadata"
params.metadata = "metadata.csv"
params.worm_model_dir = "${params.data_dir}/worm_models"

/*
~ ~ ~ > * LOG AND HELP MESSAGE SETUP
*/

if (!params.help) {
log.info '''
C E L L P R O F I L E R - N F   P I P E L I N E
===============================================
'''
    log.info ""
    log.info "Project           = ${params.project}"
    log.info "CP pipeline       = ${params.pipeline}"
    log.info "Groups            = ${params.groups}"
    log.info "Output            = ${params.out}"
    log.info "mask              = ${params.mask}"
    log.info ""
    } else {
log.info '''
C E L L P R O F I L E R - N F   P I P E L I N E
===============================================
'''
    log.info "Usage:"
    log.info "The typical command for running the pipeline is as follows:"
    log.info "nextflow run main.nf --pipeline <CellProfiler pipeline to use> --project <path to your project directory>"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--project                      The path to your project directory"
    log.info "--pipeline                     The CP pipeline to use: toxin, dauer"
    log.info ""
    log.info "Optional arguments:"
    log.info "--groups                       comma separated metadata groupings for CellProfiler, default is plate,well"
    log.info "--outdir                       Output directory to place files, default is project/Analysis-{current date}"
    log.info "--mask                         The .png file to use as a well mask in input_data/well_masks/, default is 20240618_well_mask.png"
    log.info "--help                         This usage statement."
        exit 1
    }

/*
~ ~ ~ > * WORKFLOW
*/

workflow {
    
    if("${params.pipeline}" == "dauer") {
    // configure inputs for CellProfiler ONLY FOR DAUER NOW NEED TO CHANGE MODELS IF OTHER
    config_cp = Channel.fromPath("${params.raw_pipe}")
        .combine(Channel.from("${params.metadata_dir}"))
        .combine(Channel.from("${params.metadata}"))
        .combine(Channel.from("${params.worm_model_dir}"))
        .combine(Channel.from(worm_model1)) // edit here
        .combine(Channel.from(worm_model2)) // edit here
        .combine(Channel.fromPath("${params.bin_dir}/config_CP_input_dauer.R"))
        .combine(Channel.from("${params.project}"))
        .combine(Channel.from("${params.well_mask}"))
        .combine(Channel.from("${params.groups}"))
        .combine(Channel.from("${params.edited_pipe}"))
        .combine(Channel.from("${params.out}")) | config_CP_input_dauer
        //.view()

     // Run CellProfiler
    groups = config_CP_input_dauer.out.groups_file
        .splitCsv(header:true, sep: "\t")
        .map { row ->
                [row.group, file("${row.pipeline}"), file("${row.output}")]
            }
        //.view()
    
    runCP(groups)
    
    // Preprocess CellProfiler output files
    proc_cp = runCP.out.cp_output
        .last() // This ensures that all items are emitted from runCP
        .combine(Channel.from("${params.out}"))
        .combine(Channel.from(model_name1)) // HARDCODE VARIABLE NOW MAKE DEPENDENT ON PROFILE
        .combine(Channel.from(model_name2)) // HARDCODE VARIABLE NOW MAKE DEPENDENT ON PROFILE
        .combine(Channel.fromPath("${params.bin_dir}/proc_CP_output_dauer.R"))
        //.view()

    proc_CP_output_dauer(proc_cp)
    }

    if("${params.pipeline}" == "toxin") {
    // configure inputs for CellProfiler ONLY FOR DAUER NOW NEED TO CHANGE MODELS IF OTHER
    config_cp = Channel.fromPath("${params.raw_pipe}")
        .combine(Channel.from("${params.metadata_dir}"))
        .combine(Channel.from("${params.metadata}"))
        .combine(Channel.from("${params.worm_model_dir}"))
        .combine(Channel.from(worm_model1)) // edit here
        .combine(Channel.from(worm_model2)) // edit here
        .combine(Channel.from(worm_model3)) // edit here
        .combine(Channel.from(worm_model4)) // edit here
        .combine(Channel.fromPath("${params.bin_dir}/config_CP_input_toxin.R"))
        .combine(Channel.from("${params.project}"))
        .combine(Channel.from("${params.mask}"))
        .combine(Channel.from("${params.well_mask}"))
        .combine(Channel.from("${params.groups}"))
        .combine(Channel.from("${params.edited_pipe}"))
        .combine(Channel.from("${params.out}")) | config_CP_input_toxin
        //.view()

    // Run CellProfiler
    groups = config_CP_input_toxin.out.groups_file
        .splitCsv(header:true, sep: "\t")
        .map { row ->
                [row.group, file("${row.pipeline}"), file("${row.output}")]
            }
        //.view()
    
    runCP(groups)
    
    // Preprocess CellProfiler output files
    proc_cp = runCP.out.cp_output
        .last() // This ensures that all items are emitted from runCP
        .combine(Channel.from("${params.out}"))
        .combine(Channel.from(model_name1)) // HARDCODE VARIABLE NOW MAKE DEPENDENT ON PROFILE
        .combine(Channel.from(model_name2)) // HARDCODE VARIABLE NOW MAKE DEPENDENT ON PROFILE
        .combine(Channel.from(model_name3)) // HARDCODE VARIABLE NOW MAKE DEPENDENT ON PROFILE
        .combine(Channel.from(model_name4)) // HARDCODE VARIABLE NOW MAKE DEPENDENT ON PROFILE
        .combine(Channel.fromPath("${params.bin_dir}/proc_CP_output_toxin.R")) | proc_CP_output_toxin
        //.view()
    }
}

/*
~ ~ ~ > * CONFIGURE FILES FOR CELLPROFILER
*/

process config_CP_input_dauer {
    publishDir "${params.out}/pipeline", mode: 'copy', pattern: "*.cppipe"
    publishDir "${params.out}/metadata", mode: 'copy', pattern: "metadata.csv"
    publishDir "${params.out}/groups", mode: 'copy', pattern: "groups.tsv"

    input:
        tuple file(raw_pipe), val(meta_dir), val(meta), val(model_dir), val(model1), val(model2),
        file(config_script), val(project), val(mask), val(group), val(edited_pipe), val(out)

    output:
        path "*.cppipe", emit: cp_pipeline_file
        path "metadata.csv", emit: metadata_file
        path "groups.tsv", emit: groups_file
        

    """
        # Configure the raw pipeline for CellProfiler
        awk '{gsub(/METADATA_DIR/,"${meta_dir}"); print}' ${raw_pipe} | \\
        awk '{gsub(/METADATA_CSV_FILE/,"${meta}"); print}' | \\
        awk '{gsub(/WORM_MODEL_DIR/,"${model_dir}"); print}' | \\
        awk '{gsub(/MODEL1_XML_FILE/,"${model1}"); print}' | \\
        awk '{gsub(/MODEL2_XML_FILE/,"${model2}"); print}' > pipeline.cppipe

        # Configure metadata and groups for CellProfiller with config_CP_input.R
        Rscript --vanilla ${config_script} ${project} ${mask} ${group} ${edited_pipe} ${out}

    """
}

process config_CP_input_toxin {
    publishDir "${params.out}/pipeline", mode: 'copy', pattern: "*.cppipe"
    publishDir "${params.out}/metadata", mode: 'copy', pattern: "metadata.csv"
    publishDir "${params.out}/groups", mode: 'copy', pattern: "groups.tsv"

    input:
        tuple file(raw_pipe), val(meta_dir), val(meta), val(model_dir), val(model1), val(model2), val(model3), val(model4),
        file(config_script), val(project), val(mask), val(mask_path), val(group), val(edited_pipe), val(out)

    output:
        path "*.cppipe", emit: cp_pipeline_file
        path "metadata.csv", emit: metadata_file
        path "groups.tsv", emit: groups_file
        

    """
        # Configure the raw pipeline for CellProfiler
        awk '{gsub(/METADATA_DIR/,"${meta_dir}"); print}' ${raw_pipe} | \\
        awk '{gsub(/METADATA_CSV_FILE/,"${meta}"); print}' | \\
        awk '{gsub(/WELL_MASK/,"${mask}"); print}' | \\
        awk '{gsub(/WORM_MODEL_DIR/,"${model_dir}"); print}' | \\
        awk '{gsub(/MODEL1_XML_FILE/,"${model1}"); print}' | \\
        awk '{gsub(/MODEL2_XML_FILE/,"${model2}"); print}' | \\
        awk '{gsub(/MODEL3_XML_FILE/,"${model3}"); print}' | \\
        awk '{gsub(/MODEL4_XML_FILE/,"${model4}"); print}' > pipeline.cppipe

        # Configure metadata and groups for CellProfiller with config_CP_input.R
        Rscript --vanilla ${config_script} ${project} ${mask_path} ${group} ${edited_pipe} ${out}

    """
}

/*
~ ~ ~ > * RUN CELLPROFILER
*/

process runCP {

    label "cellpro"

    input:
        tuple val(group), file(pipeline), file(output)

    output:
        stdout emit: cp_output //tuple file("*.csv"), file("*.png"), emit: cp_output

    """
        # Run cellprofiler headless
        cellprofiler -c -r -p ${pipeline} \
        -g ${group} \
        -o ${output}

    """
}

/*
~ ~ ~ > * PROCESS CELLPROFILER OUTPUTS
*/

process proc_CP_output_dauer {

    //publishDir "${params.out}/processed_data", mode: 'copy', pattern: "*.RData"
    
    input:
        tuple val(cp_output), val(out_dir), val(model_name1), val(model_name2), file(proc_CP_out_script)

    output:
        //path "*.RData", emit: cp_out_dat

    """
        # remove exisitng directories if present and make fresh
        if [ -d ${out_dir}/processed_data ]; then rm -Rf ${out_dir}/processed_data; fi
        mkdir ${out_dir}/processed_data
        if [ -d ${out_dir}/processed_images ]; then rm -Rf ${out_dir}/processed_images; fi
        mkdir ${out_dir}/processed_images
        # find .csv files, concatenate them, and write new file
        find ${out_dir}/CP_output -type f -name '${model_name1}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name1}.csv
        
        # find .csv files, concatenate them, and write new file
        find ${out_dir}/CP_output -type f -name '${model_name2}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name2}.csv
        
        # move all the output images to process_images directory END WITH /?
        find ${out_dir}/CP_output -name '*.png' -exec mv {} ${out_dir}/processed_images \\;
        # Process the CellProfiler output with proc_CP_output.R
        Rscript --vanilla ${proc_CP_out_script} ${out_dir}

    """
}

process proc_CP_output_toxin {

    //publishDir "${params.out}/processed_data", mode: 'copy', pattern: "*.RData"
    //publishDir "${params.out}/processed_data", mode: 'copy'
    
    input:
        tuple val(cp_output), val(out_dir), val(model_name1), val(model_name2),
        val(model_name3), val(model_name4), file(proc_CP_out_script)

    output:
        //path "*.RData", emit: cp_out_dat

    """
        # remove exisitng directories if present and make fresh
        if [ -d ${out_dir}/processed_data ]; then rm -Rf ${out_dir}/processed_data; fi
        mkdir ${out_dir}/processed_data
        if [ -d ${out_dir}/processed_images ]; then rm -Rf ${out_dir}/processed_images; fi
        mkdir ${out_dir}/processed_images
        
        # find .csv files, concatenate them, and write new file
        find ${out_dir}/CP_output -type f -name '${model_name1}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name1}.csv
        find ${out_dir}/CP_output -type f -name '${model_name2}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name2}.csv
        find ${out_dir}/CP_output -type f -name '${model_name3}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name3}.csv
        find ${out_dir}/CP_output -type f -name '${model_name4}.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/${model_name4}.csv
        find ${out_dir}/CP_output -type f -name 'WormObjects.csv' -print0 | xargs -0 awk 'FNR>1 || NR==1 {print}' > ${out_dir}/processed_data/WormObjects.csv
        
        # move all the output images to process_images directory
        find ${out_dir}/CP_output -name '*.png' -exec mv {} ${out_dir}/processed_images \\;
        
        # Process the CellProfiler output with proc_CP_output.R
        Rscript --vanilla ${proc_CP_out_script} ${out_dir}

    """
}

/*
~ ~ ~ > * GENERATE REPORT
*/
workflow.onComplete {

    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
    { Parameters }
    ---------------------------
    Project                                 = ${params.project}
    Pipeline Used                           = ${params.pipeline}
    Result Directory                        = ${params.out}
    """

    println summary

}
