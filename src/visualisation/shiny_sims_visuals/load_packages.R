########################
## Loading R packages ##
########################


# Function will
# - Check if packages are already installed
# - If they are not, it will install the package and its dependencies
# - Load all packages 

load_packages <- Vectorize(function(package) {
  
  # Check if package is installed, suppress warning message (only logical T/F)
  condition <- suppressWarnings(!require(package, character.only = T))
  
  # Install if not installed yet
  if (condition)
    install.packages(package, dep = T)
  
  # Load
  require(package, character.only = T)
})


# Use:
# packages <- c("tidyverse", "Cairo", "pROC")
# load_packages(packages)


