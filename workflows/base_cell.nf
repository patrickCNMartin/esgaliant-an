nextflow.enable.dsl=2

include { isMac; isLinux; getContainerImage } from "${baseDir}/lib/utils.nf"

workflow baseCellWorkflow {
    take:
    containers_ready  // Channel to ensure containers are built first
    
    main:
    // Collect all containers to ensure all are built before proceeding
    // This creates a dependency so baseCellProcess waits for containers
    baseCellProcess(containers_ready.collect())
    
    emit:
    zarr_store = baseCellProcess.out.zarr_store
    plot = baseCellProcess.out.plot
}

process baseCellProcess {
    label 'base_cell'
    
    container "${getContainerImage('base_cell')}"
    
    input:
    val containers_ready  // Dummy input to ensure containers are built first
    
    output:
    // Output paths - Nextflow will capture these from the mounted baseDir
    path "${params.zarrStore}", emit: zarr_store, type: 'dir'
    path "${params.baseCellOutput}", emit: plot, optional: true
    
    script:
    // Use absolute paths to ensure zarr_cache is accessible outside container
    // These paths are in baseDir which Nextflow automatically mounts for containers
    // The zarr store will be written inside the container but accessible outside
    // because baseDir is mounted and the output is declared above
    def zarrCache = file(params.zarrCache).toAbsolutePath().toString()
    def zarrStore = file(params.zarrStore).toAbsolutePath().toString()
    def plotOutput = params.baseCellOutput ? file(params.baseCellOutput).toAbsolutePath().toString() : null
    def command = params.baseCellCommand
    def atlas = params.baseCellAtlas
    def organism = params.baseCellOrganism
    def keepAll = params.baseCellKeepAll ? "--keep-all" : ""
    def maxChunkSize = params.baseCellMaxChunkSize
    def minChunkSize = params.baseCellMinChunkSize
    def distanceMetrics = params.baseCellDistanceMetrics ? "--distance-metrics ${params.baseCellDistanceMetrics.join(' ')}" : ""
    // Only pass --output for plot command
    def outputPath = (command == "plot" && plotOutput) ? "--output ${plotOutput}" : ""
    def maxGenes = params.baseCellMaxGenes ? "--max-genes ${params.baseCellMaxGenes}" : ""
    def maxCellTypes = params.baseCellMaxCellTypes ? "--max-cell-types ${params.baseCellMaxCellTypes}" : ""
    
    """
    # Create zarr_cache directory if it doesn't exist (using absolute path)
    mkdir -p ${zarrCache}
    
    # Run base_cell.py
    python /app/bin/base_cell.py \\
        ${command} \\
        --zarr-store ${zarrStore} \\
        --atlas ${atlas} \\
        --organism ${organism} \\
        ${keepAll} \\
        --max-chunk-size ${maxChunkSize} \\
        --min-chunk-size ${minChunkSize} \\
        ${distanceMetrics} \\
        ${outputPath} \\
        ${maxGenes} \\
        ${maxCellTypes}
    """
}

