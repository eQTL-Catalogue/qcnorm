#!/usr/bin/env nextflow
nextflow.enable.dsl=2

qtl_inputs_def_file = file("$baseDir/assets/qtlmap_inputs.tsv")
qtl_inputs_def_file.copyTo("${params.outdir}/qtlmap_inputs.tsv")
qtlmap_inputs_file = file("${params.outdir}/qtlmap_inputs.tsv")
pheno_metadata_list = [
    "ge": params.ge_pheno_meta_path,
    "exon": params.exon_pheno_meta_path,
    "tx": params.tx_pheno_meta_path,
    "txrev": params.txrev_pheno_meta_path,
    "microarray": params.array_pheno_meta_path
]

include { normalise_microarray; normalise_RNAseq_ge ; normalise_RNAseq_exon ; normalise_RNAseq_tx ; normalise_RNAseq_txrev; normalise_RNAseq_leafcutter } from  '../modules/normalisation_multi'

def add_to_qtlmap_input_tsv(qtlgroup_quantiletpm_ch, quant_method){
    // item[0] : run_id
    // item[1] : study_name
    // item[2] : sample_meta_path
    // item[3] : vcf_file
    // item[4] : median_tpm_file (for microarray normalised_splitted_count_matrix)
    // item[5] : normalised_splitted_count_matrix
    
    if (quant_method=="microarray"){
        qtlgroup_quantiletpm_ch
        .collectFile(storeDir:"${params.outdir}/qtl_group_inputs") { item ->
                    [ "${item[0]}_${item[1]}_${quant_method}_tsv_inputs.txt", 
                    "${item[1]}_${quant_method}_${item[4].baseName - ".tsv" - item[1] - '.'}\t" + 
                    "${params.outdir}/${item[0]}/${item[1]}/normalised/${quant_method}/qtl_group_split_norm/${item[4].fileName}\t" + 
                    "${pheno_metadata_list[quant_method]}\t" + 
                    "${item[2]}\t" + 
                    "${item[3]}\t" + 
                    "null.txt" + '\n' ]
                }
                .subscribe{ qtlmap_inputs_file.append(it.text) }
    }
    else if (quant_method=="leafcutter"){
        qtlgroup_quantiletpm_ch
        .collectFile(storeDir:"${params.outdir}/qtl_group_inputs") { item ->
                    [ "${item[0]}_${item[1]}_${quant_method}_tsv_inputs.txt", 
                    "${item[1]}_${quant_method}_${item[5].baseName - ".tsv" - item[1] - '.'}\t" + 
                    "${params.outdir}/${item[0]}/${item[1]}/normalised/${quant_method}/qtl_group_split_norm/${item[5].fileName}\t" + 
                    "${params.outdir}/${item[0]}/${item[1]}/normalised/${quant_method}/leafcutter_metadata.txt.gz\t" + 
                    "${item[2]}\t" + 
                    "${item[3]}\t" + 
                    "${params.outdir}/${item[0]}/${item[1]}/normalised/${item[4].fileName}" + '\n' ]
                }
                .subscribe{ qtlmap_inputs_file.append(it.text) }
    } else {
        qtlgroup_quantiletpm_ch
        .collectFile(storeDir:"${params.outdir}/qtl_group_inputs") { item ->
                    [ "${item[0]}_${item[1]}_${quant_method}_tsv_inputs.txt", 
                    "${item[1]}_${quant_method}_${item[5].baseName - ".tsv" - item[1] - '.'}\t" + 
                    "${params.outdir}/${item[0]}/${item[1]}/normalised/${quant_method}/qtl_group_split_norm/${item[5].fileName}\t" + 
                    "${pheno_metadata_list[quant_method]}\t" + 
                    "${item[2]}\t" + 
                    "${item[3]}\t" + 
                    "${params.outdir}/${item[0]}/${item[1]}/normalised/${item[4].fileName}" + '\n' ]
                }
                .subscribe{ qtlmap_inputs_file.append(it.text) }
    }
}

workflow {
    normalise()
}

workflow normalise {
        //run_id	study_name	quant_results_path	sample_meta_path	vcf_file
    Channel.fromPath(params.input_tsv)
        .ifEmpty { error "Cannot find input_tsv file in: ${params.input_tsv}" }
        .splitCsv(header: true, sep: '\t', strip: true)
        .map{row -> [ row.run_id, row.study_name, file(row.quant_results_path), file(row.sample_meta_path) ]}
        .set { study_file_ge_ch }

    Channel.fromPath(params.input_tsv)
        .ifEmpty { error "Cannot find input_tsv file in: ${params.input_tsv}" }
        .splitCsv(header: true, sep: '\t', strip: true)
        .map{row -> [ row.run_id, row.study_name, row.sample_meta_path, row.vcf_file ]}
        .set { output_tsv_ch }
        
    if (params.is_microarray){
        normalise_microarray(study_file_ge_ch, Channel.fromPath(params.array_pheno_meta_path, checkIfExists: true).collect())

        add_to_qtlmap_input_tsv(output_tsv_ch
            .combine(normalise_microarray.out.qtlmap_tsv_input_ch, by: 0).transpose(), "microarray")
    }
    else {
        normalise_RNAseq_ge(study_file_ge_ch, Channel.fromPath(params.ge_pheno_meta_path, checkIfExists: true).collect())
        
        add_to_qtlmap_input_tsv(output_tsv_ch
            .join(normalise_RNAseq_ge.out.median_tpm_file)
            .combine(normalise_RNAseq_ge.out.qtlmap_tsv_input_ch, by: 0).transpose(), "ge")
        
        if (!params.skip_exon_norm) {
            normalise_RNAseq_exon(normalise_RNAseq_ge.out.inputs_with_quant_tpm_ch, Channel.fromPath(params.exon_pheno_meta_path, checkIfExists: true).collect())
        
            add_to_qtlmap_input_tsv(output_tsv_ch
                .join(normalise_RNAseq_ge.out.median_tpm_file)
                .combine(normalise_RNAseq_exon.out.qtlmap_tsv_input_ch, by: 0).transpose(), "exon")
        }
        if (!params.skip_tx_norm) {
            normalise_RNAseq_tx(normalise_RNAseq_ge.out.inputs_with_quant_tpm_ch, Channel.fromPath(params.tx_pheno_meta_path, checkIfExists: true).collect())

            add_to_qtlmap_input_tsv(output_tsv_ch
                .join(normalise_RNAseq_ge.out.median_tpm_file)
                .combine(normalise_RNAseq_tx.out.qtlmap_tsv_input_ch, by: 0).transpose(), "tx")
        } 
        if (!params.skip_txrev_norm) {
            normalise_RNAseq_txrev(normalise_RNAseq_ge.out.inputs_with_quant_tpm_ch, Channel.fromPath(params.txrev_pheno_meta_path, checkIfExists: true).collect())

            add_to_qtlmap_input_tsv(output_tsv_ch
                .join(normalise_RNAseq_ge.out.median_tpm_file)
                .combine(normalise_RNAseq_txrev.out.qtlmap_tsv_input_ch, by: 0).transpose(), "txrev")
        }
        if (!params.skip_leafcutter_norm) {
            normalise_RNAseq_leafcutter(normalise_RNAseq_ge.out.inputs_with_quant_tpm_ch, 
                                        Channel.fromPath(params.leafcutter_transcript_meta, checkIfExists: true).collect(),
                                        Channel.fromPath(params.leafcutter_intron_annotation, checkIfExists: true).collect())


            add_to_qtlmap_input_tsv(output_tsv_ch
                .join(normalise_RNAseq_ge.out.median_tpm_file)
                .combine(normalise_RNAseq_leafcutter.out.qtlmap_tsv_input_ch, by: 0).transpose(), "leafcutter")
        }
    }
}
