#!/usr/bin/env nextflow
nextflow.preview.dsl=2

params.num_pc = 3
params.data_name = "Study_name_temp"

Channel
    .fromPath(params.vcf_file)
    .ifEmpty { exit 1, "Samples genotype vcf not found: ${params.vcf_file}" } 
    .set { vcf_file_ch }

Channel
    .fromPath(params.ref_genome)
    .ifEmpty { exit 1, "Reference genotype vcf not found: ${params.ref_genome}" } 
    .set { ref_genome_ch }

Channel
    .fromPath(params.populations_file)
    .ifEmpty { exit 1, "Populations metadata file not found: ${params.populations_file}" } 
    .set { populations_file_ch }

// if(params.exclude_population){
//     Channel
//         .fromPath(params.ids_to_remove_file)
//         .ifEmpty { exit 1, "Populations metadata file not found: ${params.ids_to_remove_file}" } 
//         .set { ids_to_remove_file_ch }
// }
    
workflow pop_assign {
    // get:
    // ref_genome_ch
    // vcf_file_ch
    // populations_file_ch

    // main:
    refVCFtoBED(ref_genome_ch)
    refBedNoDubl = params.exclude_population ? 
        removeFamilyFromRef(refVCFtoBED.out, ids_to_remove_file_ch) : 
        removeDublFromRef(refVCFtoBED.out)

    getSNPsFromRef(refBedNoDubl)
     //getSNPsFromRef.out.ref_snps_list

    sampleVCftoBED(vcf_file_ch) 

    calculate_relatedness_matrix(sampleVCftoBED.out.sample_bed, sampleVCftoBED.out.relatedness_sample_list)
    // calculate_relatedness_matrix.out
    
    // extractSharedSNPsFromSampleGen(sampleVCftoBED.out.sample_bed, 
    //     getSNPsFromRef.out.ref_snps_list)
    // // extractSharedSNPsFromSampleGen.out.sample_gen_overlapped
    // // extractSharedSNPsFromSampleGen.out.sample_gen_overlapped_snplist
    
    // extractSharedSNPsFromRefGen(refBedNoDubl,
    //     extractSharedSNPsFromSampleGen.out.sample_gen_overlapped_snplist)
    // // extractSharedSNPsFromRefGen.out.ref_overlapped
    
    // calcKinsMatrices(extractSharedSNPsFromRefGen.out.ref_overlapped)
    // // calcKinsMatrices.out.ref_overlapped_kins
    
    // calcRefPcaAndLoads(calcKinsMatrices.out.ref_overlapped_kins,
    //     extractSharedSNPsFromRefGen.out.ref_overlapped)
    // //calcRefPcaAndLoads.out.ref_overlapped_loads
    // //calcRefPcaAndLoads.out.ref_overlapped_pca

    // mapSampleGenToRef(extractSharedSNPsFromSampleGen.out.sample_gen_overlapped,
    //     calcRefPcaAndLoads.out.ref_overlapped_loads)
    // // mapSampleGenToRef.out.sample_gen_scores
    // plot_pca(mapSampleGenToRef.out.sample_gen_scores,
    //          calcRefPcaAndLoads.out.ref_overlapped_pca,
    //          populations_file_ch)
}

process refVCFtoBED{
    storeDir "$baseDir/reference_vcf_cache/"
    
    input:
    path "ref.vcf.gz"

    output:
    tuple file('ref.bed'), file('ref.bim'), file('ref.fam')

    script:
    """
    # do ld pruning and save resulted snps in file  
    plink2 --vcf ref.vcf.gz --vcf-half-call h --indep-pairwise 50000 200 0.05 --out ref_pruned_varaints_list --threads ${task.cpus}

    # make bfiles for pruned 1000 genome proj 
    plink2 --vcf ref.vcf.gz --vcf-half-call h --extract ref_pruned_varaints_list.prune.in --make-bed --out ref
    """
}

process removeFamilyFromRef{
    storeDir "$baseDir/reference_vcf_cache/"

    input:
    tuple file('ref.bed'), file('ref.bim'), file('ref.fam')
    path 'ids_to_remove.txt' 

    output:
    tuple file('ref_no_dubl.bed'), file('ref_no_dubl.bim'), file('ref_no_dubl.fam')

    script:
    """
    plink2 --bfile ref --remove-fam ids_to_remove.txt --make-bed --out ref

    # finds dublicate vars
    plink2 --bfile ref --list-duplicate-vars --out dubl

    # delete dublicate vars
    plink2 --bfile ref --exclude dubl.dupvar --snps-only --make-bed --out ref_no_dubl
    """
}

process removeDublFromRef{
    storeDir "$baseDir/reference_vcf_cache/"

    input:
    tuple file('ref.bed'), file('ref.bim'), file('ref.fam')
  
    output:
    tuple file('ref_no_dubl.bed'), file('ref_no_dubl.bim'), file('ref_no_dubl.fam')

    script:
    """
    # finds dublicate vars
    plink2 --bfile ref --list-duplicate-vars --out dubl
    
    # delete dublicate vars
    plink2 --bfile ref --exclude dubl.dupvar --snps-only --make-bed --out ref_no_dubl
    """
}

process getSNPsFromRef{
    storeDir "$baseDir/reference_vcf_cache/"

    input: 
    tuple file('ref_no_dubl.bed'), file('ref_no_dubl.bim'), file('ref_no_dubl.fam')

    output:
    path 'ref_snps_list.snplist', emit: ref_snps_list

    script:
    """
    plink2 --bfile ref_no_dubl --write-snplist --out ref_snps_list --snps-only
    """
}

// convert vcf file to plink binary file (.bed)
process sampleVCftoBED{
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    file 'sample.vcf.gz'

    output:
    tuple file('sample_genotype.bed'), file('sample_genotype.bim'), file('sample_genotype.fam'), emit: sample_bed
    path 'sample_list.txt', emit: relatedness_sample_list

    script:
    """
    #Convert VCF to plink
    plink2 --vcf sample.vcf.gz --out sample_genotype --threads ${task.cpus}
    plink2 --bfile sample_genotype --list-duplicate-vars --out list_dubl
    plink2 --bfile sample_genotype --exclude list_dubl.dupvar --snps-only --make-bed --out sample_genotype

    #Extract sample ids
    bcftools query -l sample.vcf.gz > sample_list.txt
    """
}

process calculate_relatedness_matrix{
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple file('sample_genotype.bed'), file('sample_genotype.bim'), file('sample_genotype.fam')
    path 'sample_list.txt'

    output:
    path "relatedness_matrix.tsv"

    script:
    """
    #Filter low MAF and low HWE snps
    plink2 --bfile sample_genotype --maf 0.05 --hwe 1e-6 --hwe-all --make-bed --out binary_filtered

    #Perform LD pruning
    plink2 --bfile binary_filtered --indep-pairwise 250 50 0.2 --out pruned_genotypes

    #Make pruned plink file 
    plink2 --bfile binary_filtered --extract pruned_genotypes.prune.in --make-bed --out binary_pruned

    #Calculate relatedness
    plink2 --make-rel square --bfile binary_pruned

    #Format relatedness matrix
    Rscript $baseDir/bin/pop_assign/format_kinship.R \\
        --kinship plink.rel \\
        --fam sample_list.txt \\
        --out relatedness_matrix.tsv
    """
}

// process extractSharedSNPsFromSampleGen{
//     input:
//     set file ('sample_genotype.bed'), file('sample_genotype.bim'), file('sample_genotype.fam')
//     file 'ref_snps_list.snplist'

//     output:
//     set file('sample_gen_overlapped.bed'), file('sample_gen_overlapped.bim'), file('sample_gen_overlapped.fam'), emit: sample_gen_overlapped
//     file 'overlapped_snps.snplist', emit: sample_gen_overlapped_snplist

//     // extract snps present in pruned data from new dataset
//     // --make-bed makes sure bfiles created!
//     script:
//     """
//     plink2 --bfile sample_genotype --extract ref_snps_list.snplist --make-bed --out sample_gen_overlapped
//     plink2 --bfile sample_gen_overlapped --write-snplist --out overlapped_snps
//     """
// }

// process extractSharedSNPsFromRefGen {
//     input:
//     set file('ref.bed'), file('ref.bim'), file('ref.fam')
//     file 'overlapped_snps.snplist'

//     output:
//     set file('ref_overlapped.bed'), file('ref_overlapped.bim'), file('ref_overlapped.fam'), emit: ref_overlapped

//     script:
//     """
//     plink2 --bfile ref --extract overlapped_snps.snplist --make-bed --out ref_overlapped
//     """
// }

// process calcKinsMatrices{
//     input:
//     set file('ref_overlapped.bed'), file('ref_overlapped.bim'), file('ref_overlapped.fam')

//     output:
//     set file("ref_overlapped_kins.grm.bin"), file("ref_overlapped_kins.grm.id"), file("ref_overlapped_kins.grm.adjust"), file("ref_overlapped_kins.grm.details"), emit: ref_overlapped_kins
    
//     script:
//     """
//     ldak --calc-kins-direct ref_overlapped_kins --bfile ref_overlapped --ignore-weights YES --power -0.25
//     """
// }

// process calcRefPcaAndLoads{
//     publishDir "${params.outdir}", mode: 'copy'

//     input:
//     set file("ref_overlapped_kins.grm.bin"), file("ref_overlapped_kins.grm.id"), file("ref_overlapped_kins.grm.adjust"), file("ref_overlapped_kins.grm.details")
//     set file('ref_overlapped.bed'), file('ref_overlapped.bim'), file('ref_overlapped.fam')

//     output:
//     file 'ref_overlapped_loads.load', emit: ref_overlapped_loads
//     file 'ref_overlapped_pca.vect', emit: ref_overlapped_pca
     
//     script: 
//     """
//     ldak --pca ref_overlapped_pca --grm ref_overlapped_kins --axes $params.num_pc
//     ldak --calc-pca-loads ref_overlapped_loads --grm ref_overlapped_kins --pcastem ref_overlapped_pca --bfile ref_overlapped
//     """
// }

// process mapSampleGenToRef{
//     publishDir "${params.outdir}", mode: 'copy'
    
//     input:
//     set file('sample_gen_overlapped.bed'), file('sample_gen_overlapped.bim'), file('sample_gen_overlapped.fam')
//     file 'ref_overlapped_loads.load'

//     output:
//     file 'sample_gen_scores.profile.adj', emit: sample_gen_scores

//     script:
//     """
//     ldak --calc-scores sample_gen_scores --bfile sample_gen_overlapped --scorefile ref_overlapped_loads.load --power 0
//     """
// }

// process plot_pca{
//     publishDir "${params.outdir}", mode: 'copy'
    
//     input:
//     file 'sample_gen_scores.profile.adj'
//     file 'ref_overlapped_pca.vect'
//     file 'samples_data.tsv'

//     output:
//     set file('main_pca.png'), file('projections_only.png'), file('projections_on_ref.png'), file('populations.tsv'), file('knn_threshold.png'), file('knn.png')

//     script:
//     """
//     Rscript $baseDir/bin/pop_assign/plot_pca.R ref_overlapped_pca.vect sample_gen_scores.profile.adj samples_data.tsv $params.data_name
//     """

// }

// workflow.onComplete { 
// 	println ( workflow.success ? "Done!" : "Oops ... something went wrong" )
// }

