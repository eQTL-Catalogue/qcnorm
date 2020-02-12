#!/usr/bin/env nextflow
nextflow.preview.dsl=2

exp_matrix_path = params.is_microarray ? 
                    params.exp_matrix_path : 
                    "${params.quant_results_path}/featureCounts/merged_gene_counts.txt"

if(!params.skip_ge_norm){
    Channel.fromPath(params.ge_exp_matrix_path ? 
            params.ge_exp_matrix_path : exp_matrix_path, 
            checkIfExists: true)
        .set { ge_count_matrix_ch }
}
if(!params.skip_exon_norm && !params.is_microarray){
    Channel.fromPath(params.exon_exp_matrix_path ? 
            params.exon_exp_matrix_path : "${params.quant_results_path}/dexseq_exon_counts/merged_exon_counts.tsv", 
            checkIfExists: true)
        .set { exon_count_matrix_ch }
}
if(!params.skip_tx_norm && !params.is_microarray){
    Channel.fromPath(params.tx_usage_matrix_path ? 
            params.tx_usage_matrix_path : 
            "${params.quant_results_path}/Salmon/merged_counts/TPM/gencode.v30.transcripts.TPM.merged.txt", 
            checkIfExists: true)
        .set { tx_usage_matrix_ch }
}
if(!params.skip_txrev_norm && !params.is_microarray){
    Channel.fromPath(params.txrev_usage_matrix_path ? 
            params.txrev_usage_matrix_path : 
            "${params.quant_results_path}/Salmon/merged_counts/TPM/", 
            checkIfExists: true, type: 'dir')
        .set { txrev_usage_matrix_ch }
}
if(!params.skip_leafcutter_norm && !params.is_microarray){
    Channel.fromPath(params.leafcutter_usage_matrix_path ? 
            params.leafcutter_usage_matrix_path : 
            "${params.quant_results_path}/leafcutter/leafcutter_perind_numers.counts.formatted.gz", 
            checkIfExists: true)
        .set { leafcutter_matrix_ch }
}

Channel
    .fromPath(params.sample_meta_path, checkIfExists: true)
    .set { sample_metadata_ch }


workflow normalise {
    if (params.is_microarray){
        normalise_microarray(
            ge_count_matrix_ch, 
            sample_metadata_ch, 
            Channel.fromPath(params.array_pheno_meta_path, checkIfExists: true))
    }
    else {
        if (!params.skip_ge_norm) {
            normalise_RNAseq_ge(
                ge_count_matrix_ch, 
                sample_metadata_ch, 
                Channel.fromPath(params.ge_pheno_meta_path, checkIfExists: true))
        }
        if (!params.skip_exon_norm) {
            normalise_RNAseq_exon(
                exon_count_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.exon_pheno_meta_path, checkIfExists: true))
        }
        if (!params.skip_tx_norm) {
            normalise_RNAseq_tx(
                tx_usage_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.tx_pheno_meta_path, checkIfExists: true))
        } 
        if (!params.skip_txrev_norm) {
            normalise_RNAseq_txrev(
                txrev_usage_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.txrev_pheno_meta_path, checkIfExists: true))
        }
        if (!params.skip_leafcutter_norm) {
            normalise_RNAseq_leafcutter(
                leafcutter_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.leafcutter_transcript_meta, checkIfExists: true),
                Channel.fromPath(params.leafcutter_intron_annotation, checkIfExists: true))
        }
    }

}

process normalise_microarray{
    publishDir "${params.outdir}/${params.study_name}/normalised/", mode: 'copy'
    
    label 'process_medium'
    container = 'kerimoff/eqtlutils:latest'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    
    output:
    path "*._norm_exprs.tsv"

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
    publishDir "${params.outdir}/${params.study_name}/normalised/", mode: 'copy'
    
    label 'process_medium'
    container = 'kerimoff/eqtlutils:latest'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    
    output:
    path "*.gene_counts_cqn_norm.tsv"

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
    publishDir "${params.outdir}/${params.study_name}/normalised/", mode: 'copy'
    
    label 'process_medium'
    container = 'kerimoff/eqtlutils:latest'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    
    output:
    path "*.exon_counts_cqn_norm.tsv"

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
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}

process normalise_RNAseq_tx{
    publishDir "${params.outdir}/${params.study_name}/normalised/", mode: 'copy'
    
    label 'process_medium'
    container = 'kerimoff/eqtlutils:latest'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    
    output:
    path "*_qnorm.tsv"

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
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}

process normalise_RNAseq_txrev{
    publishDir "${params.outdir}/${params.study_name}/normalised/", mode: 'copy'
    
    label 'process_medium'
    container = 'kerimoff/eqtlutils:latest'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    
    output:
    path "*_qnorm.tsv"

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
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}

process normalise_RNAseq_leafcutter{
    publishDir "${params.outdir}/${params.study_name}/normalised/", mode: 'copy'
    
    label 'process_medium'
    container = 'kauralasoo/eqtlutils:96d357d24e1b14e312298bdbd2deb0fd408660a3' //change it to more stable image
    
    input:
    path count_matrix
    path sample_metadata
    path transcript_meta
    path intron_annotation
    
    output:
    path "*_qnorm.tsv"
    path "leafcutter_metadata.txt.gz"

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
      $filter_qc\
      $keep_XY\
      $eqtl_utils_path

    """
}