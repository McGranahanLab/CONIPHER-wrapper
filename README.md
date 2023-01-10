# CONIPHER 

## CONIPHER clustering and tree building wrapper

This is a README detailing how to run both mutation clustering and phylogenetic tree building using CONIPHER (COrrecting Noise In PHylogenetic Evaluation and Reconstruction).

--- 
### Setup

To create the conda enrivonment needed to run CONIPHER clustering and tree building with one wrapper script, follow the steps below. Our R package for CONIPHER tree building is available for download [here](https://github.com/McGranahanLab/CONIPHER). Please refer to our manuscript for further details of the method.

1. Create a conda environment with correct libaries installed
```
conda create -n conipher -c bioconda -c conda-forge pyclone r-base=3.6.1 r-essentials r-tidyverse r-cowplot r-ggpubr r-fst r-biocmanager r-devtools r-seqminer
```

2. Once this has been run, activate the conda environment and start R

```
conda activate conipher
R
```

3. Subsequently install the below list of packages from an R session. NOTE: please do not update related packages during installation when prompted to do so. 

```
# Packages required for CONIPHER clustering

install.packages("mclust")
BiocManager::install("GenomicRanges")
BiocManager::install("Rsamtools")
install.packages("gplots")
install.packages("gdata")
install.packages("future")
install.packages("optparse")
install.packages("bootstrap")
#install.packages("TeachingDemos")
#devtools::install_version("NORMT3", version = "1.0.4")
#devtools::install_github("genome/bmm")
BiocManager::install("copynumber")
devtools::install_version("sequenza", version = "2.1.2")
install.packages("coin")
install.packages("wordcloud")


# CONIPHER treebuilding R package
devtools::install_github("McGranahanLab/CONIPHER")
```

4. Once all of these have been installed quit R and deactivate the conda environment

```
q()
conda deactivate
```

--- 

### Quickstart
#### Running CONIPHER clustering + tree building

We provide a wrapper bash script to run CONIPHER clustering and tree building end-to-end. To run this from the conda environment set up as above on the example case CRUKTOY001 provided, run the following command in your terminal:

```
sh 0_runningClusteringTreeBuilding.sh
```
