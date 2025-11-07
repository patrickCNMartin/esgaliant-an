nextflow.enable.dsl=2
//=============================================================================
// IMPORTS
//=============================================================================
include { isMac; isLinux; getContainerImage } from './lib/utils.nf'
include { checkAndBuildContainers } from './workflows/build_containers.nf'
include { baseCellWorkflow } from './workflows/base_cell.nf'
//=============================================================================
// WORKFLOW
//=============================================================================
workflow {
    println "Check container images and build if abscent" 
    checkAndBuildContainers()
    
    println "Running base cell workflow"
    // Pass container outputs to ensure containers are built before base cell workflow starts
    baseCellWorkflow(checkAndBuildContainers.out.containers)
    
    baseCellWorkflow.out.zarr_store.view { "Zarr store: ${it}" }
    baseCellWorkflow.out.plot.view { "Plot: ${it}" }
}