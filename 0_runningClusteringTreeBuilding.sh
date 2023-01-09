#!/bin/bash

################################################################################## Input parameters
###################################################################################################

case_id="LTXTOY003"
scriptDir=`pwd`"/"
inputTSV=${scriptDir}"/toy_tsv.tsv"
outDir=${scriptDir}"/"${case_id}


############################################################### Running clustering and treebuilding
###################################################################################################

source activate pyclone_2

scDir=${outDir}"/Clustering/"
treeDir=${outDir}"/TreeBuilding/"

mkdir -p ${scDir}
mkdir -p ${treeDir}

Rscript ${scriptDir}run_clustering.R --case_id ${case_id} --script_dir ${scriptDir} --input_tsv ${inputTSV} --working_dir ${scDir} --pyclone_version sc_nosciclone --nProcs 2
Rscript ${scriptDir}run_treebuilding.R --input_tsv ${scDir}${case_id}".SCoutput.CLEAN.tsv" --out_dir ${treeDir} --script_dir ${scriptDir} --prefix LTX

# conda deactivate


############################################################################################### End
###################################################################################################
