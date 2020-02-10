#!/usr/bin/env nextflow

nextflow.preview.dsl=2
include './modules/pop_assign' params(params)
include './modules/quality_control' params(params)

workflow {
    // pop_assign()
    println params
    quality_control()
}