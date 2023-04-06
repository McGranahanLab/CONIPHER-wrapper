# CONIPHER 

## CONIPHER clustering and tree building wrapper

This is a README detailing how to run both mutation clustering and phylogenetic tree building using CONIPHER (COrrecting Noise In PHylogenetic Evaluation and Reconstruction). NOTE: our R package for CONIPHER tree building is available for download [here](https://github.com/McGranahanLab/CONIPHER). For full details of all the inputs and expected outputs for CONIPHER clustering and tree building, refer to our protocol (https://doi.org/10.21203/rs.3.pex-2158/v1).

--- 
### Setup

Clone the github repo using the following command from your terminal and enter the directory:
```
git clone https://github.com/McGranahanLab/CONIPHER-wrapper/
cd CONIPHER-wrapper
```

#### Create CONIPHER conda environment
To be able to run CONIPHER clustering and tree building with one wrapper script, follow the steps below. 

To create the conda environment to successfully install and run CONIPHER clustering and tree building, please manually build the conda environment using the instructions below.

On your terminal, ensure you are located in the `CONIPHER-wrapper` directory.

1. Create the conda environment with the following libaries installed
```
conda create -n conipher -c conda-forge -c bioconda pyclone r-base=3.6.1 r-essentials r-tidyverse r-cowplot r-ggpubr r-fst r-biocmanager r-devtools r-seqminer
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

We provide a wrapper bash script to run CONIPHER clustering and tree building end-to-end. To run this from the conda environment set up as above on the example case CRUKTOY001 provided, first ensure you are in the `CONIPHER-wrapper` folder on your terminal, then enter the conda environment and run the clustering+treebuilding wrapper bash script as follows:

```
source activate conipher
sh 0_runningClusteringTreeBuilding.sh
```


#### Running CONIPHER tree building

We additionally provide a wrapper script to run CONIPHER tree building by itself. To run this from the conda environment set up as above on the example case CRUKTOY001 provided, first ensure you are in the `CONIPHER-wrapper` folder on your terminal, then run the following command:

```
sh 0_runningTreeBuilding.sh
```

