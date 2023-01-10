#!/bin/bash

################################################################################## Input parameters
###################################################################################################

case_id="CRUKTOY001"
scriptDir=`pwd`"/src/"
inputTSV="data/input_tsv.tsv"
outDir="data/results"


############################################################### Running clustering and treebuilding
###################################################################################################

source activate conipher

clusteringDir=${outDir}"/Clustering/"
treeDir=${outDir}"/TreeBuilding/"

mkdir -p ${clusteringDir}
mkdir -p ${treeDir}

Rscript ${scriptDir}run_clustering.R \
--case_id ${case_id} \
--script_dir ${scriptDir} \
--input_tsv ${inputTSV} \
--working_dir ${clusteringDir} \
--pyclone_version sc_nosciclone \
--nProcs 2

Rscript ${scriptDir}run_treebuilding.R \
--input_tsv ${clusteringDir}${case_id}".SCoutput.CLEAN.tsv" \
--out_dir ${treeDir} \
--script_dir ${scriptDir} \
--prefix CRUK

# conda deactivate


############################################################################################### End
###################################################################################################
