nextflow.enable.dsl=2
include { isMac; isLinux; getContainerImage } from "${baseDir}/lib/utils.nf"

workflow checkAndBuildContainers {
    
    Channel
        .fromPath("${params.containerDir}/*", type: 'dir')
        .set { container_dirs }
    
    buildContainers(container_dirs)
}

process buildContainers {
    
    input:
    each path(containerDir)
    
    output:
    path "*.{sif,tar}", optional: true
    
    script:
    def name = containerDir.name
    def cacheDir = params.containerCache
    def useMac = isMac()
    """
    mkdir -p ${cacheDir}
    
    if [ "${useMac}" = "true" ]; then
        imagePath="${cacheDir}/${name}.tar"
        
        if [ ! -f "\${imagePath}" ]; then
            echo "Building Docker image for ${name}..."
            docker build -t ${name}:latest -f ${containerDir}/Dockerfile ${containerDir}
            docker save ${name}:latest -o \${imagePath}
            echo "Built Docker image for ${name}"
        else
            echo "Docker image for ${name} already exists"
        fi
        ln -s \${imagePath} ${name}.tar
    else
        imagePath="${cacheDir}/${name}.sif"
        
        if [ ! -f "\${imagePath}" ]; then
            echo "Building Singularity image for ${name}..."
            singularity build \${imagePath} ${containerDir}/Singularity.def
            echo "Built Singularity image for ${name}"
        else
            echo "Singularity image for ${name} already exists"
        fi
        ln -s \${imagePath} ${name}.sif
    fi
    """
}