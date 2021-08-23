#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process normalise_microarray{
    publishDir "${params.publishDir}/microarray", mode: 'copy'
    
    label 'process_medium'
    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    
    output:
    path "*_norm_exprs.tsv"
    path "qtl_group_split_norm/*", emit: qtlmap_tsv_input_ch

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $count_matrix\
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
    publishDir "${params.publishDir}/ge", mode: 'copy',
        saveAs: {filename -> filename.indexOf("_tpm.tsv.gz") > 0 ? "../$filename" : "$filename"}

    label 'process_medium'
    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    
    output:
    path "*.gene_counts_cqn_int_norm.tsv"
    path "*_95quantile_tpm.tsv.gz", emit: quantile_tpm_file
    path "*_median_tpm.tsv.gz", emit: median_tpm_file
    path "qtl_group_split_norm/*", emit: qtlmap_tsv_input_ch

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $count_matrix\
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
    publishDir "${params.publishDir}/exon", mode: 'copy'
    
    label 'process_high'
    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    path tpm_quantile
    
    output:
    path "*.exon_counts_cqn_int_norm.tsv"
    path "qtl_group_split_norm/*", emit: qtlmap_tsv_input_ch

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $count_matrix\
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
    publishDir "${params.publishDir}/tx", mode: 'copy'
    
    label 'process_medium'
    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    path tpm_quantile
    
    output:
    path "*_qnorm.tsv"
    path "qtl_group_split_norm/*", emit: qtlmap_tsv_input_ch

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $count_matrix\
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
    publishDir "${params.publishDir}/txrev", mode: 'copy'
    
    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    path tpm_quantile
    
    output:
    path "*_qnorm.tsv"
    path "qtl_group_split_norm/*", emit: qtlmap_tsv_input_ch

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $count_matrix\
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

process normalise_RNAseq_leafcutter{
    publishDir "${params.publishDir}/leafcutter", mode: 'copy'
    
    label 'process_medium'
    container = 'quay.io/eqtlcatalogue/eqtlutils:v20.04.1'
    
    input:
    path count_matrix
    path sample_metadata
    path transcript_meta
    path intron_annotation
    path tpm_quantile
    
    output:
    path "*_qnorm.tsv"
    path "leafcutter_metadata.txt.gz"
    path "qtl_group_split_norm/*", emit: qtlmap_tsv_input_ch
    // TODO: Test leafcutter behaviour. never tested

    script:
    filter_qc = params.norm_filter_qc ? "--filter_qc TRUE" : ""
    keep_XY = params.norm_keep_XY ? "--keep_XY TRUE" : ""
    eqtl_utils_path = params.eqtl_utils_path ? "--eqtlutils ${params.eqtl_utils_path}" : ""
    """
    # Make leafcutter phenotype metadata file
    Rscript $baseDir/bin/normalisation/makeLeafcutterMetadata.R\
      -c $count_matrix\
      -t $transcript_meta\
      -i $intron_annotation\
      -o leafcutter_metadata.txt.gz\
      $eqtl_utils_path

    Rscript $baseDir/bin/normalisation/normaliseCountMatrix.R\
      -c $count_matrix\
      -s $sample_metadata\
      -p leafcutter_metadata.txt.gz\
      -o .\
      -q leafcutter\
      -t $tpm_quantile\
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}
