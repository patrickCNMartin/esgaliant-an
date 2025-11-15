nextflow.enable.dsl=2

include { isMac; isLinux; getContainerImage } from "${baseDir}/lib/utils.nf"

workflow baseCellWorkflow {
    take:
    containers_ready  // Channel to ensure containers are built first
    
    main:
    baseCellProcess(containers_ready.collect())
    
    emit:
    zarr_store = baseCellProcess.out.zarr_store
    plot = baseCellProcess.out.plot
}

process baseCellProcess {
    label 'base_cell'
    publishDir params.baseCellOutput, copy:true, overwrite:true
    container "${getContainerImage('base_cell')}"
    
    input:
    val containers_ready
    
    output:
    path "${params.zarrStore}", emit: zarr_store, type: 'dir'
    path "${params.plot_out}/*.png", emit: plots
    
    script:
    """
    python /app/bin/base_cell.py \
        --zarr-store "${params.zarrStore}" \
        --atlas ${params.baseCellAtlas} \
        --organism ${params.baseCellOrganism} \
        ${params.baseCellKeepAll ? '--keep-all' : ''} \
        --max-chunk-size ${params.baseCellMaxChunkSize} \
        --min-chunk-size ${params.baseCellMinChunkSize} \
        --distance-metrics ${params.baseCellDistanceMetrics} \
        --output ${params.plot_out} \
        --max-genes ${params.baseCellMaxGenes} \
        --max-cell-types ${params.baseCellMaxCellTypes}
    """
}