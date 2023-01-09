#!/bin/bash

################################################################################## Input parameters
###################################################################################################

patient="LTXTOY003"
scriptDir=`pwd`"/"
inputTSV=${scriptDir}"/toy_tsv.tsv"
outDir=${scriptDir}"/"${patient}


############################################################### Running clustering and treebuilding
###################################################################################################
source activate pyclone_2

scDir=${outDir}"/Clustering/"
treeDir=${outDir}"/TreeBuilding/"

mkdir -p ${scDir}
mkdir -p ${treeDir}

Rscript ${scriptDir}run_treebuilding.R --input_tsv ${inputTSV} --out_dir ${treeDir} --script_dir ${scriptDir} --prefix LTX

# conda deactivate


############################################################################################### End
###################################################################################################
