#!/usr/bin/env nextflow

nextflow.preview.dsl=2
include './modules/pop_assign' params(params)

workflow {
    pop_assign()
}