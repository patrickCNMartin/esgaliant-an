//=============================================================================
// HELPER FUNCTION
//=============================================================================
def isMac() { System.properties['os.name'].toLowerCase().contains('mac') }

def isLinux() { System.properties['os.name'].toLowerCase().contains('linux') }

def getContainerImage(name) {
    def imagePath = isMac ? "${params.containerCache}/${name}.tar" 
                          : "${params.containerCache}/${name}.sif"
    if (!file(imagePath).exists()) {
        error "Container image not found: ${imagePath}. Run buildContainers first."
    }
    return imagePath
}
