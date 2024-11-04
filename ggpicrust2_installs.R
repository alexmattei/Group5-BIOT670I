### This script contains the download commands for all of the libraries needed.
### Run it before loading the libraries for other scripts.
### Commands are taken from the ggpicrust2 github page: https://github.com/cafferychen777/ggpicrust2


install.packages("ggpicrust2")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("phyloseq", "ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "BiocManager", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# Install ggfun, needed for the pca plots
install.packages('ggfun', repos = c('https://yulab-smu.r-universe.dev', 'https://cloud.r-project.org'))

