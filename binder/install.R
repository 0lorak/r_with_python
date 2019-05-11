### install regular packages

install.packages("reticulate") # python support in RMarkdown
install.packages("ggplot2") # for plotting
install.packages(c("rmarkdown", "caTools", "bitops")) # for knitting
install.packages("lattice") # to be able to install markovchain package
install.packages("Matrix") # to be able to install markovchain package
install.packages("RcppParallel") # to be able to install markovchain package
install.packages("matlab") # to be able to install markovchain package
install.packages("markovchain") # to be able to work with Markov Chains

### install bioconductor packages
# install.packages("BiocManager")
# BiocManager::install("package")

### install GitHub packages (tag = commit, branch or release tag)
# install.packages("devtools")
# devtools::install_github("user/repo", ref = "tag")
