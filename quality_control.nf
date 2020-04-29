#!/usr/bin/env nextflow
nextflow.preview.dsl=2

exp_matrix_path = params.is_microarray ? params.microarray_exp_matrix_path : 
                    params.quant_results_path ? "${params.quant_results_path}/featureCounts/merged_gene_counts.txt" : params.ge_exp_matrix_path

if (!params.is_microarray) {
    params.mbv_path = params.quant_results_path ? "${params.quant_results_path}/MBV" : 
                        params.mbv_files_dir ? params.mbv_files_dir : { exit 1, "Error: Please provide --mbv_files_dir or --quant_results_path parameter. "}
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

if(params.pop_assign_projections){
    Channel
        .fromPath(params.pop_assign_projections, checkIfExists: true)
        .set { pop_assign_projections_custom_ch }
}

workflow {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 
    if(params.pop_assign_projections=="") {exit 1, "Error: Please provide --pop_assign_projections parameter. "} 

    build_qc_report(
        script_path_ch,
        ge_count_matrix_ch,
        sample_metadata_ch,
        pheno_metadata_ch,
        pop_assign_projections_custom_ch
    )
}

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

process build_qc_report{
    publishDir "${params.outdir}/${params.study_name}/QC", mode: 'copy'
    
    container = 'eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    path script_path
    path ge_count_matrix
    path sample_metadata
    path pheno_metadata
    path pop_assign_projections
    
    output:
    path "*.html"
    path "*.rds"
    path "*matrix.tsv"

    script:
    mbv_path = params.is_microarray ? "" : "mbv_files_dir = \"${params.mbv_path}\","
    eqtl_utils_path = params.eqtl_utils_path ? "eqtl_utils_path = \"${params.eqtl_utils_path}\"," : ""
    """
    #!/usr/bin/env Rscript

    work_dir = getwd()

    rmarkdown::render(
    input = "$script_path", 
    output_dir = work_dir,
    output_file = "${params.study_name}_QC_report.html",
    params = list(
        work_dir = work_dir,
        set_title = "${params.study_name} QC report", 
        set_author = "${params.author_name}",
        count_matrix_path = "$ge_count_matrix",
        sample_meta_path = "$sample_metadata",
        phenotype_meta_path = "$pheno_metadata", $mbv_path $eqtl_utils_path
        projections = "$pop_assign_projections"))
    """
    
}
