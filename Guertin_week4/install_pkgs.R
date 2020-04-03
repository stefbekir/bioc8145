# installing packages on Rivanna 
# Set the path to the Parallel libraries
local_lib <- "~/R/parallel_lib/3.6"
if (! dir.exists(local_lib)){
   dir.create(local_lib, recursive=TRUE)
}
.libPaths('~/R/parallel_lib/3.6')

# Set up CRAN repository for all packages
fav_repo <- 'https://mirrors.nics.utk.edu/cran'
local({
   r <- getOption("repos")
   r['CRAN'] <- fav_repo
   options(repos=r)
})

# Start list of CRAN packages that need to be installed
pkgs <- c( 'digest', 
           'rlang',
           'Rcpp',
           'htmltools',
           'tidyselect',
           'glue',
           'ellipsis',
           'fansi',
           'vctrs',
           'tibble',
          'dplyr', 
          'BiocManager',
          'pkgload'
	 )
# Loop through list of CRAN packages
# If not installed already, install it
installed <- rownames(installed.packages())
for (p in pkgs){
  if ( ! p %in% installed){
    install.packages(p)
  }
}

# Install bioconductor packages
BiocManager::install("S4Vectors")
BiocManager::install("IRanges")
BiocManager::install("DelayedArray")
BiocManager::install("SummarizedExperiment")
BiocManager::install("IRanges")
BiocManager::install("DESeq2")
BiocManager::install("fgsea")
BiocManager::install("gage")
BiocManager::install("limma")

install.packages('lattice')
install.packages('dplyr')
install.packages('ggplot2')

#*********************************************************************#
