nextflow.enable.dsl=2
include { isMac; isLinux; getContainerImage } from "${baseDir}/lib/utils.nf"

workflow checkAndBuildContainers {
    
    Channel
        .fromPath("${params.containerDir}/*", type: 'dir')
        .set { container_dirs }
    
    buildContainers(container_dirs)
    
    emit:
    containers = buildContainers.out
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
            # Change to containerDir so Docker can find pyproject.toml in build context
            cd ${containerDir}
            docker build -t ${name}:latest -f Dockerfile .
            docker save ${name}:latest -o \${imagePath}
            echo "Built Docker image for ${name}"
        else
            echo "Docker image for ${name} already exists"
        fi
        ln -s \${imagePath} ${name}.tar
    else
        imagePath="${cacheDir}/${name}.sif"
        
        if [ ! -f "\${imagePath}" ]; then
            echo "Building apptainer image for ${name}..."
            module load apptainer
            # Change to containerDir so apptainer can find pyproject.toml
            cd ${containerDir}
            apptainer build --fakeroot \${imagePath} ${name}.def
            echo "Built apptainer image for ${name}"
        else
            echo "Apptainer image for ${name} already exists"
        fi
        ln -s \${imagePath} ${name}.sif
    fi
    """
}