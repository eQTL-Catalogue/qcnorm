#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.publishDir="${params.outdir}/${params.study_name}/normalised"
qtl_inputs_def_file = file("$baseDir/assets/qtlmap_inputs.tsv")
qtl_inputs_def_file.copyTo("${params.outdir}/${params.study_name}/${params.study_name}_qtlmap_inputs.tsv")
qtlmap_inputs_file = file("${params.outdir}/${params.study_name}/${params.study_name}_qtlmap_inputs.tsv")
pheno_metadata_list = [
    "ge": params.ge_pheno_meta_path,
    "exon": params.exon_pheno_meta_path,
    "tx": params.tx_pheno_meta_path,
    "txrev": params.txrev_pheno_meta_path
]

//group_id	quant_results_path	sample_meta_path	outdir
Channel.fromPath(params.input_tsv)
    .ifEmpty { error "Cannot find input_tsv file in: ${params.input_tsv}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.study_name, file(row.quant_results_path), file(row.sample_meta_path), file(row.vcf_file), row.outdir]}
    .view()
    .set { study_file_ge_ch }

include { normalise_RNAseq_ge ; normalise_RNAseq_exon ; normalise_RNAseq_tx ; normalise_RNAseq_txrev } from  '../modules/normalisation_multi'

// def add_to_qtlmap_input_tsv(qtlgroup_quantiletpm_ch, quant_method){
//     if (quant_method=="microarray"){
//         qtlgroup_quantiletpm_ch
//         .collectFile(storeDir:"${params.outdir}/${params.study_name}/qtl_group_inputs") { item ->
//                     [ "${params.study_name}_${quant_method}_tsv_inputs.txt", 
//                     "${params.study_name}_${quant_method}_${item.baseName - params.study_name - '.'}\t" + 
//                     "${params.publishDir}/${quant_method}/qtl_group_split_norm/${item.fileName}\t" + 
//                     "${params.array_pheno_meta_path}\t" + 
//                     "${params.sample_meta_path}\t" + 
//                     "${params.vcf_file}\t" + 
//                     "null.txt" + '\n' ]
//                 }
//                 .subscribe{ qtlmap_inputs_file.append(it.text) }
//     }
//     else {
//         qtlgroup_quantiletpm_ch
//         .collectFile(storeDir:"${params.outdir}/${params.study_name}/qtl_group_inputs") { item ->
//                     [ "${params.study_name}_${quant_method}_tsv_inputs.txt", 
//                     "${params.study_name}_${quant_method}_${item[0].baseName - params.study_name - '.'}\t" + 
//                     "${params.publishDir}/${quant_method}/qtl_group_split_norm/${item[0].fileName}\t" + 
//                     "${pheno_metadata_list[quant_method]}\t" + 
//                     "${params.sample_meta_path}\t" + 
//                     "${params.vcf_file}\t" + 
//                     "${params.publishDir}/${item[1].fileName}" + '\n' ]
//                 }
//                 .subscribe{ qtlmap_inputs_file.append(it.text) }
//     }
// }

workflow {
    normalise()
}

workflow normalise {
    normalise_RNAseq_ge(study_file_ge_ch, Channel.fromPath(params.ge_pheno_meta_path, checkIfExists: true).collect())
    // add_to_qtlmap_input_tsv(normalise_RNAseq_ge.out.qtlmap_tsv_input_ch.flatten()
    //     .combine(normalise_RNAseq_ge.out.median_tpm_file), "ge")
    if (!params.skip_exon_norm) {
        normalise_RNAseq_exon(normalise_RNAseq_ge.out.inputs_with_quant_tpm_ch, Channel.fromPath(params.exon_pheno_meta_path, checkIfExists: true).collect())

        // add_to_qtlmap_input_tsv(normalise_RNAseq_exon.out.qtlmap_tsv_input_ch.flatten()
        //     .combine(median_tpm_file_ch), "exon")
    }
    if (!params.skip_tx_norm) {
        normalise_RNAseq_tx(normalise_RNAseq_ge.out.inputs_with_quant_tpm_ch, Channel.fromPath(params.tx_pheno_meta_path, checkIfExists: true).collect())

        // add_to_qtlmap_input_tsv(normalise_RNAseq_tx.out.qtlmap_tsv_input_ch.flatten()
        //     .combine(median_tpm_file_ch), "tx")
    } 
    if (!params.skip_txrev_norm) {
        normalise_RNAseq_txrev(normalise_RNAseq_ge.out.inputs_with_quant_tpm_ch, Channel.fromPath(params.txrev_pheno_meta_path, checkIfExists: true).collect())

        // add_to_qtlmap_input_tsv(normalise_RNAseq_txrev.out.qtlmap_tsv_input_ch.flatten()
        //     .combine(median_tpm_file_ch), "txrev")
    }
}
