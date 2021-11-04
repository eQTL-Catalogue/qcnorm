#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process normalise_microarray{
    publishDir "${params.outdir}/$run_id/$study_name/normalised/microarray", mode: 'copy', pattern: "qtl_group_split_norm/*"
    
    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    tuple val(run_id), val(study_name), file(quant_results_path), file(sample_metadata)
    path pheno_metadata
    
    output:
    tuple val(run_id), file("qtl_group_split_norm/*"), emit: qtlmap_tsv_input_ch

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $quant_results_path\
      -s $sample_metadata\
      -p $pheno_metadata\
      -o .\
      -q HumanHT-12_V4\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}

process normalise_RNAseq_ge{
    publishDir "${params.outdir}/$run_id/$study_name/normalised/ge", mode: 'copy', pattern: "norm_not_filtered/*"
    publishDir "${params.outdir}/$run_id/$study_name/normalised/ge", mode: 'copy', pattern: "qtl_group_split_norm/*"
    publishDir "${params.outdir}/$run_id/$study_name/normalised/", mode: 'copy', pattern: "*_tpm.tsv.gz"

    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    tuple val(run_id), val(study_name), file(quant_results_path), file(sample_metadata)
    path pheno_metadata
    
    output:
    path "norm_not_filtered/*"
    path "*_95quantile_tpm.tsv.gz", emit: quantile_tpm_file
    tuple val(run_id), file("*_median_tpm.tsv.gz"), emit: median_tpm_file
    tuple val(run_id), file("qtl_group_split_norm/*"), emit: qtlmap_tsv_input_ch
    tuple val(run_id), val(study_name), file(quant_results_path), file(sample_metadata), file("*_95quantile_tpm.tsv.gz"), emit: inputs_with_quant_tpm_ch

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $quant_results_path/featureCounts/merged_gene_counts.txt\
      -s $sample_metadata\
      -p $pheno_metadata\
      -o .\
      -q gene_counts\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}

process normalise_RNAseq_exon{
    publishDir "${params.outdir}/$run_id/$study_name/normalised/exon", mode: 'copy'
    
    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    tuple val(run_id), val(study_name), file(quant_results_path), file(sample_metadata), file(tpm_quantile)
    path pheno_metadata

    output:
    path "norm_not_filtered/*"
    tuple val(run_id), file("qtl_group_split_norm/*"), emit: qtlmap_tsv_input_ch

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $quant_results_path/dexseq_exon_counts/merged_exon_counts.tsv\
      -s $sample_metadata\
      -p $pheno_metadata\
      -o .\
      -q exon_counts\
      -t $tpm_quantile\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}

process normalise_RNAseq_tx{
    publishDir "${params.outdir}/$run_id/$study_name/normalised/tx", mode: 'copy'
    
    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    tuple val(run_id), val(study_name), file(quant_results_path), file(sample_metadata), file(tpm_quantile)
    path pheno_metadata
    
    output:
    path "norm_not_filtered/*"
    tuple val(run_id), file("qtl_group_split_norm/*"), emit: qtlmap_tsv_input_ch

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $quant_results_path/Salmon/merged_counts/TPM/gencode.v30.transcripts.TPM.merged.txt\
      -s $sample_metadata\
      -p $pheno_metadata\
      -o .\
      -q transcript_usage\
      -t $tpm_quantile\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}

process normalise_RNAseq_txrev{
    publishDir "${params.outdir}/$run_id/$study_name/normalised/txrev", mode: 'copy'
    
    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    tuple val(run_id), val(study_name), file(quant_results_path), file(sample_metadata), file(tpm_quantile)
    path pheno_metadata
    
    output:
    path "norm_not_filtered/*"
    tuple val(run_id), file("qtl_group_split_norm/*"), emit: qtlmap_tsv_input_ch

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $quant_results_path/Salmon/merged_counts/TPM/\
      -s $sample_metadata\
      -p $pheno_metadata\
      -o .\
      -q txrevise\
      -t $tpm_quantile\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}
