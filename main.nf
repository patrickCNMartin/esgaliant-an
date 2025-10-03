nextflow.enable.dsl=2
//=============================================================================
// IMPORTS
//=============================================================================
include { isMac; isLinux; getContainerImage } from './lib/utils.nf'
include { checkAndBuildContainers } from './workflows/build_containers.nf'
//=============================================================================
// WORKFLOW
//=============================================================================
workflow {
    println "Check container images and build if abscent" 
    checkAndBuildContainers()
}