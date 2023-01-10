#!/bin/bash

################################################################################## Input parameters
###################################################################################################

case_id="CRUKTOY001"
scriptDir=`pwd`"/src/"
inputTSV=`pwd`"/data/input_tsv.tsv"
outDir=`pwd`"/data/results"

###################################################################### Running treebuilding wrapper
###################################################################################################

source activate conipher

treeDir=${outDir}"/TreeBuilding/"

mkdir -p ${treeDir}

Rscript ${scriptDir}run_treebuilding.R \
--input_tsv ${inputTSV} \
--out_dir ${treeDir} \
--script_dir ${scriptDir} \
--prefix CRUK


############################################################################################### End
###################################################################################################
