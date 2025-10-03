nextflow.enable.dsl=2
workflow buildContainers {

    process buildContainers {

        output:
        path "${params.containerCache}/*", emit: built_images

        script:
        """
        mkdir -p ${params.containerCache}
        
        for container in ${params.containerDir}/*; do
            name=\$(basename \$container)
            
            if [ "${isMac}" = "true" ]; then
                defFile="\$container/Dockerfile"
                imagePath="${params.containerCache}/\${name}.tar"
                if [ ! -f "\$imagePath" ]; then
                    docker build -t \${name}:latest -f \$defFile \$container
                    docker save \${name}:latest -o \$imagePath
                    echo "Built Docker image for \$name"
                else
                    echo "Docker image for \$name already exists"
                fi
            else
                defFile="\$container/Singularity.def"
                imagePath="${params.containerCache}/\${name}.sif"
                if [ ! -f "\$imagePath" ]; then
                    singularity build \$imagePath \$defFile
                    echo "Built Singularity image for \$name"
                else
                    echo "Singularity image for \$name already exists"
                fi
            fi
        done
        """
    }

    emit:
    buildContainers.built_images
}
