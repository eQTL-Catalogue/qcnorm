#!/usr/bin/env nextflow
nextflow.enable.dsl=2

exp_matrix_path = params.is_microarray ? params.microarray_exp_matrix_path : 
                    params.quant_results_path ? "${params.quant_results_path}/featureCounts/merged_gene_counts.txt" : params.ge_exp_matrix_path

if (!params.is_microarray) {
    params.mbv_path = params.quant_results_path ? "${params.quant_results_path}/MBV" : params.mbv_files_dir 
}
Channel.fromPath(exp_matrix_path, checkIfExists: true)
        .set { ge_count_matrix_ch }

Channel
    .fromPath(params.sample_meta_path, checkIfExists: true)
    .set { sample_metadata_ch }

pheno_meta_path = params.is_microarray ? params.array_pheno_meta_path : params.ge_pheno_meta_path
Channel
    .fromPath(pheno_meta_path, checkIfExists: true)
    .set { pheno_metadata_ch }

script_path = params.is_microarray ? "$baseDir/bin/QC/array_study_QC_report.Rmd" : "$baseDir/bin/QC/rnaseq_study_QC_report.Rmd"
Channel
    .fromPath(script_path, checkIfExists: true)
    .set { script_path_ch }

include { build_qc_report } from  '../modules/quality_control'

// workflow {
//     if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 
//     if(params.pop_assign_projections=="") {exit 1, "Error: Please provide --pop_assign_projections parameter. "} 

//     build_qc_report(
//         script_path_ch,
//         ge_count_matrix_ch,
//         sample_metadata_ch,
//         pheno_metadata_ch,
//         pop_assign_projections_custom_ch
//     )
// }

workflow quality_control {
    take: 
    pop_assign_projections_ch
    
    main:
    build_qc_report(
        script_path_ch,
        ge_count_matrix_ch,
        sample_metadata_ch,
        pheno_metadata_ch,
        pop_assign_projections_ch
    )
}
