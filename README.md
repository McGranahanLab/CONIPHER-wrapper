# CONIPHER 

## CONIPHER clustering and tree building wrapper

This is a README detailing how to run both mutation clustering and phylogenetic tree building using CONIPHER (COrrecting Noise In PHylogenetic Evaluation and Reconstruction). NOTE: our R package for CONIPHER tree building is available for download [here](https://github.com/McGranahanLab/CONIPHER). Please refer to our manuscript (XXX) for further details of the method.

--- 
### Setup

Clone the github repo using the following command from your terminal and enter the directory:
```
git clone https://github.com/McGranahanLab/CONIPHER-wrapper/
cd CONIPHER-wrapper
```

#### Create CONIPHER conda environment
To be able to run CONIPHER clustering and tree building with one wrapper script, follow the steps below. 

To create the conda environment to successfully install and run CONIPHER clustering and tree building, either import the CONIPHER conda environment from the .yaml file provided (preferred) or manually build the conda environment using the instructions below.

##### VERSION 1 - create conda environment from yaml file

On your terminal, ensure you are located in the `CONIPHER-wrapper` directory and create the conda environment by entering the following command:
```
conda env create -f conipher_env.yml
```


##### VERSION 2 - create conda environment manually

1. Create the conda environment with the following libaries installed
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

We provide a wrapper bash script to run CONIPHER clustering and tree building end-to-end. To run this from the conda environment set up as above on the example case CRUKTOY001 provided, first ensure you are in the `CONIPHER-wrapper` folder on your terminal, then run the following command:

```
sh 0_runningClusteringTreeBuilding.sh
```


#### Running CONIPHER tree building

We additionally provide a wrapper script to run CONIPHER tree building by itself. To run this from the conda environment set up as above on the example case CRUKTOY001 provided, first ensure you are in the `CONIPHER-wrapper` folder on your terminal, then run the following command:

```
sh 0_runningTreeBuilding.sh
```

