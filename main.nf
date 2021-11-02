#!/usr/bin/env nextflow
nextflow.enable.dsl=2

workflow pop_assign_only {
    include {pop_assign} from './workflows/pop_assign_wf' params(params)
    
    pop_assign()
}

workflow qc_only {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 
    
    include {quality_control} from './workflows/quality_control_wf' params(params)

    quality_control(Channel.fromPath(params.pop_assign_projections, checkIfExists: true))
}

workflow norm_only {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 

    include {normalise} from './workflows/normalisation_wf' params(params)

    normalise()
}

workflow pop_assign_and_qc {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 
    
    include {pop_assign} from './workflows/pop_assign_wf' params(params)
    include {quality_control} from './workflows/quality_control_wf' params(params)

    pop_assign()
    quality_control(pop_assign.out.projections_ch)
}

workflow qc_and_norm {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 

    include {quality_control} from './workflows/quality_control_wf' params(params)
    include {normalise} from './workflows/normalisation_wf' params(params)

    quality_control(Channel.fromPath(params.pop_assign_projections, checkIfExists: true))
    normalise()
}

workflow norm_only_with_tsv {
    include {normalise} from './workflows/normalisation_multi_wf' params(params)

    normalise()
}

workflow {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 

    include {pop_assign} from './workflows/pop_assign_wf' params(params)
    include {quality_control} from './workflows/quality_control_wf' params(params)
    include {normalise} from './workflows/normalisation_wf' params(params)
    
    pop_assign()
    quality_control(pop_assign.out.projections_ch)
    normalise()
}
