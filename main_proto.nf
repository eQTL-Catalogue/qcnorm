#!/usr/bin/env nextflow

nextflow.preview.dsl=2

// include './modules/pop_assign' params(params)
// include './modules/quality_control' params(params)
include './modules/normalisation' params(params)

workflow {
    if(params.study_name=="") {exit 1, "Error: Please provide --study_name parameter. "} 
    
    // pop_assign()
    // quality_control()
    normalise()
}