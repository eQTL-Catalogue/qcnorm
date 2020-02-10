#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.data_name = "Study_name_temp"
params.pop_assign_projections = "${params.outdir}/pop_assign/projections_comb.tsv"

Channel
    .fromPath("${params.quant_results_path}/featureCounts/merged_gene_counts.txt", checkIfExists: true)
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

Channel
    .fromPath(params.pop_assign_projections, checkIfExists: true)
    .set { pop_assign_projections_ch }


workflow quality_control {
    build_qc_report(
        script_path_ch,
        ge_count_matrix_ch,
        sample_metadata_ch,
        pheno_metadata_ch,
        pop_assign_projections_ch
    )
}

process build_qc_report{
    publishDir "${params.outdir}/QC", mode: 'copy'
    
    label 'process_low'
    container = 'kerimoff/qc_report:latest'
    
    input:
    path script_path
    path ge_count_matrix
    path sample_metadata
    path pheno_metadata
    path pop_assign_projections
    
    output:
    path "*.html"

    script:
    mbv_path = params.is_microarray ? "" : "mbv_files_dir = \"${params.quant_results_path}/MBV\","
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
