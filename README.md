# CONIPHER 

## CONIPHER clustering and tree building wrapper

This is a README detailing how to run both mutation clustering and phylogenetic tree building using CONIPHER (COrrecting Noise In PHylogenetic Evaluation and Reconstruction). NOTE: our R package for CONIPHER is available for download [here](https://github.com/McGranahanLab/CONIPHER). For full details of all the inputs and expected outputs for CONIPHER clustering and tree building, refer to our protocol (https://doi.org/10.21203/rs.3.pex-2158/v1).

--- 
### Setup

Clone the github repo using the following command from your terminal and enter the directory:
```
git clone git@github.com:McGranahanLab/CONIPHER-wrapper.git
cd CONIPHER-wrapper
```

#### Create CONIPHER conda environment
To be able to run CONIPHER clustering and tree building with one wrapper script, follow the steps below. 

To create the conda environment to successfully install and run CONIPHER clustering and tree building, please manually build the conda environment using the instructions below.

On your terminal, ensure you are located in the `CONIPHER-wrapper` directory.

1. Create the conda environment with the following libaries installed
```
conda create -n conipher -c conda-forge -c bioconda conipher
```

2. Once this has been run, activate the conda environment

```
conda activate conipher
```

--- 

### Quickstart
#### Running CONIPHER clustering + tree building

We provide a wrapper bash script to run CONIPHER clustering and tree building end-to-end. To run this from the conda environment set up as above on the example case CRUKTOY001 provided, first ensure you are in the `CONIPHER-wrapper` folder on your terminal, then enter the conda environment and run the clustering+treebuilding wrapper bash script as follows:

```
sh wrapper_conipher.sh
```

If running CONIPHER wrapper on one individual tumour case, specify the case identifier, desired output directory, location of the input .tsv file within the ``wrapper_conipher.sh`` script.



#### Running CONIPHER clustering and tree building separately

We additionally provide a wrapper script to run CONIPHER clustering and CONIPHER tree building stages individually. To run these from the conda environment set up as above on the example case CRUKTOY001 provided, first ensure you are in the `CONIPHER-wrapper` folder on your terminal, then run the following command for clustering:

```
sh wrapper_clustering.sh
```

Run the following command for tree building:

```
sh wrapper_treebuilding.sh
```

Similarly to above, the parameters can be changed within the bash scripts to run CONIPHER stages on individual cases.
