#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.publishDir="${params.outdir}/${params.study_name}/normalised"
qtl_inputs_def_file = file("$baseDir/assets/qtlmap_inputs.tsv")
qtl_inputs_def_file.copyTo("${params.outdir}/qtlmap_inputs.tsv")
qtlmap_inputs_file = file("${params.outdir}/qtlmap_inputs.tsv")
pheno_metadata_list = [
    "ge": params.ge_pheno_meta_path,
    "exon": params.exon_pheno_meta_path,
    "tx": params.tx_pheno_meta_path,
    "txrev": params.txrev_pheno_meta_path
]

//run_id	study_name	quant_results_path	sample_meta_path	vcf_file	outdir
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

include { normalise_RNAseq_ge ; normalise_RNAseq_exon ; normalise_RNAseq_tx ; normalise_RNAseq_txrev } from  '../modules/normalisation_multi'

// qtl_subset	count_matrix	pheno_meta	sample_meta	vcf	tpm_file	covariates_file
// GTEx_ge_adipose_subcutaneous	/gpfs/space/projects/eQTLCatalogue/qcnorm/GTEx_v8_norm/group_1/GTEx/normalised/ge/qtl_group_split_norm/GTEx.adipose_subcutaneous.tsv	/gpfs/space/projects/genomic_references/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz	/gpfs/space/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/GTEx.tsv	/gpfs/space/projects/GTEx/genotypes/processed/GTEx.DS_added.filtered.vcf.gz	/gpfs/space/projects/eQTLCatalogue/qcnorm/GTEx_v8_norm/group_1/GTEx/normalised/GTEx_median_tpm.tsv.gz
def add_to_qtlmap_input_tsv(qtlgroup_quantiletpm_ch, quant_method){
    // item[0] : run_id
    // item[1] : study_name
    // item[2] : sample_meta_path
    // item[3] : vcf_file
    // item[4] : median_tpm_file
    // item[5] : normalised_splitted_count_matrix
    
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
        .collectFile(storeDir:"${params.outdir}/qtl_group_inputs") { item ->
                    [ "${item[0]}_${item[1]}_${quant_method}_tsv_inputs.txt", 
                    "${item[1]}_${quant_method}_${item[5].baseName - item[1] - '.'}\t" + 
                    "${params.outdir}/${item[0]}/${item[1]}/normalised/${quant_method}/qtl_group_split_norm/${item[5].fileName}\t" + 
                    "${pheno_metadata_list[quant_method]}\t" + 
                    "${item[2]}\t" + 
                    "${item[3]}\t" + 
                    "${params.outdir}/${item[0]}/${item[1]}/normalised/${item[5].fileName}" + '\n' ]
                }
                .subscribe{ qtlmap_inputs_file.append(it.text) }
    }
}

workflow {
    normalise()
}

workflow normalise {
    normalise_RNAseq_ge(study_file_ge_ch, Channel.fromPath(params.ge_pheno_meta_path, checkIfExists: true).collect())
    
    add_to_qtlmap_input_tsv(output_tsv_ch
        .join(normalise_RNAseq_ge.out.median_tpm_file)
        .combine(normalise_RNAseq_ge.out.qtlmap_tsv_input_ch.flatten()).view(), "ge")
    
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
