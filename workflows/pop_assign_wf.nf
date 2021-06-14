#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.num_pc = params.num_pc == null ? 3 : params.num_pc

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

if(params.exclude_population){
    Channel
        .fromPath(params.ids_to_remove_file)
        .ifEmpty { exit 1, "Populations metadata file not found: ${params.ids_to_remove_file}" } 
        .set { ids_to_remove_file_ch }
}

include { refVCFtoBED; removeFamilyFromRef; removeDublFromRef; getSNPsFromRef; sampleVCftoBED; calculateRelatednessMatrix; extractSharedSNPsFromSampleGen; extractSharedSNPsFromRefGen; calcKinsMatrices; calcRefPcaAndLoads; mapSampleGenToRef; plotPCA} from '../modules/pop_assign'

workflow {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 
    
    pop_assign()
}
    
workflow pop_assign {
    main:
        refVCFtoBED(ref_genome_ch)
        refBedNoDubl = params.exclude_population ? 
            removeFamilyFromRef(refVCFtoBED.out, ids_to_remove_file_ch) : 
            removeDublFromRef(refVCFtoBED.out)

        getSNPsFromRef(refBedNoDubl)
        // getSNPsFromRef.out.ref_snps_list

        sampleVCftoBED(vcf_file_ch)
        // sampleVCftoBED.out.sample_bed, 
        // sampleVCftoBED.out.relatedness_sample_list 

        calculateRelatednessMatrix(sampleVCftoBED.out.sample_bed, sampleVCftoBED.out.relatedness_sample_list)
        // calculateRelatednessMatrix.out
        
        extractSharedSNPsFromSampleGen(sampleVCftoBED.out.sample_bed, 
            getSNPsFromRef.out.ref_snps_list)
        // extractSharedSNPsFromSampleGen.out.sample_gen_overlapped
        // extractSharedSNPsFromSampleGen.out.sample_gen_overlapped_snplist
        
        extractSharedSNPsFromRefGen(refBedNoDubl,
            extractSharedSNPsFromSampleGen.out.sample_gen_overlapped_snplist)
        // extractSharedSNPsFromRefGen.out.ref_overlapped
        
        calcKinsMatrices(extractSharedSNPsFromRefGen.out.ref_overlapped)
        // calcKinsMatrices.out.ref_overlapped_kins
        
        calcRefPcaAndLoads(calcKinsMatrices.out.ref_overlapped_kins,
            extractSharedSNPsFromRefGen.out.ref_overlapped)
        // calcRefPcaAndLoads.out.ref_overlapped_loads
        // calcRefPcaAndLoads.out.ref_overlapped_pca

        mapSampleGenToRef(extractSharedSNPsFromSampleGen.out.sample_gen_overlapped,
            calcRefPcaAndLoads.out.ref_overlapped_loads)
        // mapSampleGenToRef.out.sample_gen_scores
        plotPCA(mapSampleGenToRef.out.sample_gen_scores,
                calcRefPcaAndLoads.out.ref_overlapped_pca,
                populations_file_ch)

    emit:
        projections_ch = plotPCA.out.pop_assign_projections_ch
}
