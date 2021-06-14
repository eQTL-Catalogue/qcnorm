#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process build_qc_report{
    publishDir "${params.outdir}/${params.study_name}/QC", mode: 'copy'
    
    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
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
