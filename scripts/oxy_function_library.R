# Reads text file with a list of packages
# Calls 'pkgTest' to check if package exists,
# downloads package if needed, then loads
loadPackages <- function(packages_to_load){
  list_tmp <- read.csv(packages_to_load)
  list_in  <- list_tmp$data.table
  for( pkg in list_in ) { 
    pkgTest(pkg) 
  }
}

pkgTest <- function(x) {
  if (!require(x, character.only = TRUE )) {
    install.packages(x,dep=TRUE)
    if (!require(x, character.only = TRUE)) {
      print("Package not found")
    } 
  } else {
    print(x)
    library(x, character.only = TRUE)
  }
}



