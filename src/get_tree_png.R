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
parser$add_argument("tree_file")
parser$add_argument("output_name")
args <- parser$parse_args()
tree_file <- args$tree_file
output_name <- args$output_name

png(output_name)
tree <- read.tree(tree_file)
plot(tree)
dev.off()
