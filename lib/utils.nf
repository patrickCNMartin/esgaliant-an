//=============================================================================
// HELPER FUNCTION
//=============================================================================
def isMac() { System.properties['os.name'].toLowerCase().contains('mac') }

def isLinux() { System.properties['os.name'].toLowerCase().contains('linux') }

def getContainerImage(name) {
    def isMacOS = isMac()
    // On Mac with Docker, use the image name (tagged during build)
    // On Linux with Singularity, use the .sif file path
    if (isMacOS) {
        // For Docker, use the image name that was tagged during build
        // The tar file is saved but Docker image is also tagged as name:latest
        return "${name}:latest"
    } else {
        // For Singularity, use the .sif file path
        def imagePath = "${params.containerCache}/${name}.sif"
        if (!file(imagePath).exists()) {
            error "Container image not found: ${imagePath}. Run buildContainers first."
        }
        return imagePath
    }
}
