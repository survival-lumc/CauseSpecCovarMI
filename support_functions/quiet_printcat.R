#########################################################
## Stopping unwanted print/cat from external functions ##
#########################################################


# Function that:
# - Takes a function/expression that by default uses print() or cat(), and stops it
# - We use it for smcfcs() and for MIcombine()
# - Source: http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html


quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 