#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {pop_assign} from './workflows/pop_assign_wf' params(params)
include {quality_control} from './workflows/quality_control_wf' params(params)
include {normalise} from './workflows/normalisation_wf' params(params)

workflow pop_assign_only {
    pop_assign()
}

workflow qc_only {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 
    
    quality_control(Channel.fromPath(params.pop_assign_projections, checkIfExists: true))
}

workflow norm_only {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 

    normalise()
}

workflow pop_assign_and_qc {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 

    pop_assign()
    quality_control(pop_assign.out.projections_ch)
}

workflow qc_and_norm {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 

    quality_control(Channel.fromPath(params.pop_assign_projections, checkIfExists: true))
    normalise()
}
workflow {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 
    
    pop_assign()
    quality_control(pop_assign.out.projections_ch)
    normalise()
}
