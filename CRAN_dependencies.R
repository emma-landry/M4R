### PACKAGES DEPENDENCIES OBTAINED BY sessionInfo()

list_of_packages <- c(
  "rstan",
  "rstanarm",
  "statmod",
  "tidyselect",
  "purrr",
  "lattice",
  "V8",
  "colorspace",
  "vctrs",
  "generics",
  "getopt",
  "loo",
  "utf8",
  "rlang",
  "pkgbuild",
  "pillar",
  "nloptr",
  "glue",
  "withr",
  "matrixStats",  
  "lifecycle",
  "plyr",
  "munsell",       
  "gtable",
  "codetools",
  "inline",
  "callr",
  "ps",
  "curl",
  "fansi",
  "bayesplot",  
  "rstantools",
  "Rcpp",
  "scales",
  "RcppParallel",
  "StanHeaders",
  "jsonlite",
  "lme4",
  "gridExtra",
  "hms",
  "processx",
  "rprojroot",
  "rbibutils",
  "Rdpack",
  "cli",
  "magrittr",
  "tibble",    
  "crayon",
  "pkgconfig",
  "MASS",     
  "ellipsis",
  "Matrix",
  "prettyunits",   
  "ggridges",
  "assertthat",
  "minqa", 
  "R6",
  "bootnlme",
  "devtools"   # Necessary for epidemia installation 
)
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]

#=== loaded via a namespace (and not attached):
install.packages(new_packages, repos = "https://cran.ma.imperial.ac.uk")

#PACKAGES CALLED DIRECTLY
list_of_packages <- c("data.table",
                      "ggplot2",
                      "optparse",
                      "readr",
                      "here",
                      "dplyr")

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
install.packages(new_packages, repos = "https://cran.ma.imperial.ac.uk")

if (!require(epidemia)){
	devtools::install_github("ImperialCollegeLondon/epidemia", ref = "master")
}


