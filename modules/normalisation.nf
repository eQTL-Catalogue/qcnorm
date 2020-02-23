#!/usr/bin/env nextflow
nextflow.preview.dsl=2
params.publishDir="${params.outdir}/${params.study_name}/normalised"
qtl_inputs_def_file = file("$baseDir/assets/qtlmap_inputs.tsv")
qtl_inputs_def_file.copyTo("${params.outdir}/${params.study_name}/${params.study_name}_qtlmap_inputs.tsv")
qtlmap_inputs_file = file("${params.outdir}/${params.study_name}/${params.study_name}_qtlmap_inputs.tsv")

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

qtlmap_tsv_input_ch = Channel.empty()

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

        normalise_RNAseq_ge.out.qtlmap_tsv_input_ch
            .combine(normalise_RNAseq_ge.out.quantile_tpm_file)
            .collectFile(storeDir:"${params.publishDir}/qtl_group_inputs") { item ->
                [ "${params.study_name}_ge_tsv_inputs.txt", 
                "${params.study_name}_${item[0].simpleName}_ge\t" + 
                "${params.publishDir}/ge/qtl_group_split_norm/${item[0].fileName}\t" + 
                "${params.ge_pheno_meta_path}\t" + 
                "${params.sample_meta_path}\t" + 
                "${params.vcf_file}\t" + 
                "${params.publishDir}/${item[1].fileName}" + '\n' ]
            }
            .subscribe{ qtlmap_inputs_file.append(it.text) }
                
        }
        tpm_quantile_ch = params.external_quantile_tpm ? 
                    Channel.fromPath(params.external_quantile_tpm, checkIfExists: true) :
                    normalise_RNAseq_ge.out.quantile_tpm_file

        if (!params.skip_exon_norm) {
            normalise_RNAseq_exon(
                exon_count_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.exon_pheno_meta_path, checkIfExists: true),
                tpm_quantile_ch)

            // qtlmap_tsv_input_ch = 
            //     Channel
            //         .from("exon", "${params.publishDir}", "${params.exon_pheno_meta_path}")
            //         .collect()
            //         .combine(normalise_RNAseq_exon.out.qtlmap_tsv_input_ch.flatten())
            //         .combine(tpm_quantile_ch.flatten())
            //         .mix(qtlmap_tsv_input_ch)
        }
        if (!params.skip_tx_norm) {
            normalise_RNAseq_tx(
                tx_usage_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.tx_pheno_meta_path, checkIfExists: true),
                tpm_quantile_ch)

            // qtlmap_tsv_input_ch = 
            //     Channel
            //         .from("tx", "${params.publishDir}", "${params.tx_pheno_meta_path}")
            //         .collect()
            //         .combine(normalise_RNAseq_tx.out.qtlmap_tsv_input_ch.flatten())
            //         .combine(tpm_quantile_ch.flatten())
            //         .mix(qtlmap_tsv_input_ch)
        } 
        if (!params.skip_txrev_norm) {
            normalise_RNAseq_txrev(
                txrev_usage_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.txrev_pheno_meta_path, checkIfExists: true),
                tpm_quantile_ch)

            // qtlmap_tsv_input_ch = 
            //     Channel
            //         .from("txrev", "${params.publishDir}", "${params.txrev_pheno_meta_path}")
            //         .collect()
            //         .combine(normalise_RNAseq_txrev.out.qtlmap_tsv_input_ch.flatten())
            //         .combine(tpm_quantile_ch.flatten())
            //         .mix(qtlmap_tsv_input_ch)
        }
        if (!params.skip_leafcutter_norm) {
            normalise_RNAseq_leafcutter(
                leafcutter_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.leafcutter_transcript_meta, checkIfExists: true),
                Channel.fromPath(params.leafcutter_intron_annotation, checkIfExists: true),
                tpm_quantile_ch)

            // qtlmap_tsv_input_ch
            //     .mix(normalise_RNAseq_leafcutter.out.qtlmap_tsv_input_ch)
            //     .set { qtlmap_tsv_input_ch }
        }
    }
    // generate_qtlmap_input_tsv(qtlmap_tsv_input_ch, Channel.fromPath("$baseDir/assets/qtlmap_inputs.tsv"))
}

process normalise_microarray{
    publishDir "${params.outdir}/${params.study_name}/normalised/miroarray", mode: 'copy'
    
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
    publishDir "${params.outdir}/${params.study_name}/normalised/ge", mode: 'copy',
        saveAs: {filename -> filename.indexOf("_tpm.tsv.gz") > 0 ? "../$filename" : "$filename"}

    label 'process_medium'
    container = 'kerimoff/eqtlutils:latest'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    
    output:
    path "*.gene_counts_cqn_norm.tsv"
    path "*_95quantile_tpm.tsv.gz", emit: quantile_tpm_file
    path "*_median_tpm.tsv.gz"
    path "qtl_group_split_norm/*", emit: qtlmap_tsv_input_ch
    // tuple val("ge"), val("${params.outdir}/${params.study_name}/normalised"), path("qtl_group_split_norm/*"), val(params.ge_pheno_meta_path), path("*_95quantile_tpm.tsv.gz"), emit: qtlmap_tsv_input_ch
    
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
    publishDir "${params.outdir}/${params.study_name}/normalised/exon", mode: 'copy'
    
    label 'process_medium'
    container = 'kerimoff/eqtlutils:latest'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    path tpm_quantile
    
    output:
    path "*.exon_counts_cqn_norm.tsv"
    path "qtl_group_split_norm/*", emit: qtlmap_tsv_input_ch
    // tuple val("exon"), val("${params.outdir}/${params.study_name}/normalised"), path("qtl_group_split_norm/*"), val(params.exon_pheno_meta_path), path(tpm_quantile), emit: qtlmap_tsv_input_ch


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
    publishDir "${params.outdir}/${params.study_name}/normalised/tx", mode: 'copy'
    
    label 'process_medium'
    container = 'kerimoff/eqtlutils:latest'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    path tpm_quantile
    
    output:
    path "*_qnorm.tsv"
    path "qtl_group_split_norm/*", emit: qtlmap_tsv_input_ch
    // tuple val("tx"), val("${params.outdir}/${params.study_name}/normalised"), path("qtl_group_split_norm/*"), val(params.tx_pheno_meta_path), path(tpm_quantile), emit: qtlmap_tsv_input_ch

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
    publishDir "${params.outdir}/${params.study_name}/normalised/txrev", mode: 'copy'
    
    label 'process_medium'
    container = 'kerimoff/eqtlutils:latest'
    
    input:
    path count_matrix
    path sample_metadata
    path pheno_metadata
    path tpm_quantile
    
    output:
    path "*_qnorm.tsv"
    path "qtl_group_split_norm/*", emit: qtlmap_tsv_input_ch
    // tuple val("txrev"), val("${params.outdir}/${params.study_name}/normalised"), path("qtl_group_split_norm/*"), val(params.txrev_pheno_meta_path), path(tpm_quantile), emit: qtlmap_tsv_input_ch

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
    publishDir "${params.outdir}/${params.study_name}/normalised/leafcutter", mode: 'copy'
    
    label 'process_medium'
    container = 'kauralasoo/eqtlutils:96d357d24e1b14e312298bdbd2deb0fd408660a3' //change it to more stable image
    
    input:
    path count_matrix
    path sample_metadata
    path transcript_meta
    path intron_annotation
    path tpm_quantile
    
    output:
    path "*_qnorm.tsv"
    path "leafcutter_metadata.txt.gz"
    path "qtl_group_split_norm/"
    // val "${params.outdir}/${params.study_name}/normalised/leafcutter/qtl_group_split_norm/*", emit: qtlmap_tsv_input_ch
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

// process generate_qtlmap_input_tsv {
//     publishDir "${params.outdir}/${params.study_name}/", mode: 'copy'
//     stageInMode "link"

//     label 'process_low'
//     container = 'kerimoff/eqtlutils:latest'
    
//     input:
//     path input_tsv_file
//     path "qtlmap_inputs.tsv"

//     output:
//     path "qtlmap_inputs.tsv"

//     script:
//     """
//     echo -e '${published_file.baseName}.$quant_method\t$published_dir/$quant_method/${published_file.fileName}\t$pheno_metadata_path\t${params.sample_meta_path}\t${params.vcf_file}\t$published_dir/${quantile_tpm_file.fileName}\n'\
//      >> qtlmap_inputs.tsv
//     """

// }