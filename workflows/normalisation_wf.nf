#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.publishDir="${params.outdir}/${params.study_name}/normalised"
qtl_inputs_def_file = file("$baseDir/assets/qtlmap_inputs.tsv")
qtl_inputs_def_file.copyTo("${params.outdir}/${params.study_name}/${params.study_name}_qtlmap_inputs.tsv")
qtlmap_inputs_file = file("${params.outdir}/${params.study_name}/${params.study_name}_qtlmap_inputs.tsv")
pheno_metadata_list = [
    "ge": params.ge_pheno_meta_path,
    "exon": params.exon_pheno_meta_path,
    "tx": params.tx_pheno_meta_path,
    "txrev": params.txrev_pheno_meta_path
]

exp_matrix_path = params.is_microarray ? 
                    params.microarray_exp_matrix_path : 
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



include { normalise_microarray ; normalise_RNAseq_ge ; normalise_RNAseq_exon ; normalise_RNAseq_tx ; normalise_RNAseq_txrev ; normalise_RNAseq_leafcutter} from  '../modules/normalisation'


def add_to_qtlmap_input_tsv(qtlgroup_quantiletpm_ch, quant_method){
    if (quant_method=="microarray"){
        qtlgroup_quantiletpm_ch
        .collectFile(storeDir:"${params.outdir}/${params.study_name}/qtl_group_inputs") { item ->
                    [ "${params.study_name}_${quant_method}_tsv_inputs.txt", 
                    "${params.study_name}_${quant_method}_${item.baseName - params.study_name - '.'}\t" + 
                    "${params.publishDir}/${quant_method}/qtl_group_split_norm/${item.fileName}\t" + 
                    "${params.array_pheno_meta_path}\t" + 
                    "${params.sample_meta_path}\t" + 
                    "${params.vcf_file}\t" + 
                    "null.txt" + '\n' ]
                }
                .subscribe{ qtlmap_inputs_file.append(it.text) }
    }
    else {
        qtlgroup_quantiletpm_ch
        .collectFile(storeDir:"${params.outdir}/${params.study_name}/qtl_group_inputs") { item ->
                    [ "${params.study_name}_${quant_method}_tsv_inputs.txt", 
                    "${params.study_name}_${quant_method}_${item[0].baseName - params.study_name - '.'}\t" + 
                    "${params.publishDir}/${quant_method}/qtl_group_split_norm/${item[0].fileName}\t" + 
                    "${pheno_metadata_list[quant_method]}\t" + 
                    "${params.sample_meta_path}\t" + 
                    "${params.vcf_file}\t" + 
                    "${params.publishDir}/${item[1].fileName}" + '\n' ]
                }
                .subscribe{ qtlmap_inputs_file.append(it.text) }
    }
}

workflow {
    normalise()
}

workflow normalise {
    if (params.is_microarray){
        normalise_microarray(
            ge_count_matrix_ch, 
            sample_metadata_ch, 
            Channel.fromPath(params.array_pheno_meta_path, checkIfExists: true))

        add_to_qtlmap_input_tsv(normalise_microarray.out.qtlmap_tsv_input_ch.flatten(), "microarray")
    }
    else {
        if (!params.skip_ge_norm) {
            normalise_RNAseq_ge(
                ge_count_matrix_ch, 
                sample_metadata_ch, 
                Channel.fromPath(params.ge_pheno_meta_path, checkIfExists: true))
        
        add_to_qtlmap_input_tsv(normalise_RNAseq_ge.out.qtlmap_tsv_input_ch.flatten()
            .combine(normalise_RNAseq_ge.out.median_tpm_file), "ge")
                
        }
        tpm_quantile_ch = params.external_quantile_tpm ? 
                    Channel.fromPath(params.external_quantile_tpm, checkIfExists: true) :
                    normalise_RNAseq_ge.out.quantile_tpm_file

        median_tpm_file_ch = params.external_median_tpm ? 
                    Channel.fromPath(params.external_median_tpm, checkIfExists: true) :
                    normalise_RNAseq_ge.out.median_tpm_file

        if (!params.skip_exon_norm) {
            normalise_RNAseq_exon(
                exon_count_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.exon_pheno_meta_path, checkIfExists: true),
                tpm_quantile_ch)

            add_to_qtlmap_input_tsv(normalise_RNAseq_exon.out.qtlmap_tsv_input_ch.flatten()
                .combine(median_tpm_file_ch), "exon")
        }
        if (!params.skip_tx_norm) {
            normalise_RNAseq_tx(
                tx_usage_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.tx_pheno_meta_path, checkIfExists: true),
                tpm_quantile_ch)

            add_to_qtlmap_input_tsv(normalise_RNAseq_tx.out.qtlmap_tsv_input_ch.flatten()
                .combine(median_tpm_file_ch), "tx")
        } 
        if (!params.skip_txrev_norm) {
            normalise_RNAseq_txrev(
                txrev_usage_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.txrev_pheno_meta_path, checkIfExists: true),
                tpm_quantile_ch)

            add_to_qtlmap_input_tsv(normalise_RNAseq_txrev.out.qtlmap_tsv_input_ch.flatten()
                .combine(median_tpm_file_ch), "txrev")
        }
        if (!params.skip_leafcutter_norm) {
            normalise_RNAseq_leafcutter(
                leafcutter_matrix_ch,
                sample_metadata_ch,
                Channel.fromPath(params.leafcutter_transcript_meta, checkIfExists: true),
                Channel.fromPath(params.leafcutter_intron_annotation, checkIfExists: true),
                tpm_quantile_ch)

            // add_to_qtlmap_input_tsv(normalise_RNAseq_leafcutter.out.qtlmap_tsv_input_ch.flatten()
            //     .combine(median_tpm_file_ch), "leafcutter")
        }
    }
}
