nextflow.enable.dsl=2
//=============================================================================
// IMPORTS
//=============================================================================
include { isMac, isLinux, getContainerImage } from './lib/utils.nf'
include { buildContainers } from './workflows/build_containers.nf'
//=============================================================================
// WORKFLOW
//=============================================================================
workflow {
    println "Hello World" 
    built_images = buildContainers()
}