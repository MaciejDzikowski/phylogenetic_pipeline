if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
required_packages <- c("ape", "argparse", "phangorn")
for (package in required_packages) { #Installs packages if not yet installed
    if (!requireNamespace(package, quietly=TRUE))
        BiocManager::install(package)
}

library(ape)
library(argparse)
library(phangorn)


parser <-ArgumentParser()
parser$add_argument("nj_trees_file")
parser$add_argument("output_dir")
args <- parser$parse_args()
nj_trees_file <- args$nj_trees_file
output_dir <- args$output_dir

treelist <- read.tree(file=nj_trees_file)  # ./data/nj_trees.nwk
consensus_tree <- consensus(treelist, p=0.5)
output_name <- paste0(output_dir, "/consensus_tree.nwk")
write.tree(consensus_tree, file=output_name)
# plot(consensus_tree)
