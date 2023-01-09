cat('\n======== RUNNING TREE BUILDING ===========\n')

#### =========== PACKAGES REQUIRED ========= ####
cat('\nLoading required packages')
require(optparse)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(CONIPHER))

#### =========== SETUP ========= ####
option_list = list(
  
  # generic options
  make_option(c("--input_tsv"), type="character", default=NULL,
              help="location of input tsv", metavar="character"),
  make_option(c("--out_dir"), type="character", default=NULL,
              help="where should the output be saved", metavar="character"),
  make_option(c("--script_dir"), type="character", default=NULL,
              help="where are the scripts that should be sourced located?", metavar="character"),

  #specific options for tree building parameters
  make_option(c("--ccf_buffer"), type="integer", default=10,
              help="buffer used for CCF calculations", metavar="character"),
  make_option(c("--prefix"), type="character", default=NULL,
              help="sample prefix", metavar="character"),
  make_option(c("--pval_cutoff"), type="double", default=0.01,
              help="p-value for copy number calculation", metavar="character"),
  make_option(c("--use_boot"), type="logical", default="TRUE",
              help="TRUE or FALSE should bootstrapping be used", metavar="character"),
  make_option(c("--merge_clusters"), type="logical", default="TRUE",
              help="TRUE or FALSE should similar clusters be merged if possible", metavar="character"),
  make_option(c("--correct_cpn_clusters"), type="logical", default="TRUE",
              help="TRUE or FALSE should clusters driven by copy number be removed", metavar="character"),
  make_option(c("--adjust_noisy_clusters"), type="logical", default="FALSE",
              help="TRUE or FALSE should noisy clusters be adjusted?", metavar="character"),
  make_option(c("--adjust_noisy_clusters_prop"), type="double", default=0.05,
              help="what is minimum proportion of mutations should be present in a region to avoid cluster adjustment", metavar="character"),
  make_option(c("--min_ccf"), type="double", default=0.01,
              help="what is minimum CCF to consider a mutation as present? [range 0-1]", metavar="character"),
  make_option(c("--min_cluster_size"), type="double", default=5,
              help="what is the minimum number of mutations in a cluster for it to be included in analysis?", metavar="character"),
  make_option(c("--multi_trees"), type="logical", default="TRUE",
              help="explore alternative trees", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input_tsv               <- opt$input_tsv
prefix                  <- opt$prefix
script_dir              <- opt$script_dir
out_dir                 <- opt$out_dir
run.multi.trees         <- opt$multi_trees

# assign parameters from opt
ccf_buffer                    <- opt$ccf_buffer 
pval_cutoff                   <- opt$pval_cutoff 
use_boot                      <- opt$use_boot 
merge_clusters                <- opt$merge_clusters 
correct_cpn_clusters          <- opt$correct_cpn_clusters 
adjust_noisy_clusters         <- opt$adjust_noisy_clusters 
adjust_noisy_clusters_prop    <- opt$adjust_noisy_clusters_prop 
min_ccf                       <- opt$min_ccf 
min_cluster_size              <- opt$min_cluster_size 

# create output directories
if (!file.exists(out_dir)) { 
  dir.create(out_dir, showWarnings = TRUE, recursive = TRUE, mode = "0775")
}

# load the tsv file
input_table <- read.table(input_tsv, sep = "\t", stringsAsFactors = FALSE, header = TRUE)

#### =========== PREOCESS INPUT DATA ========= ####

# preprocess input data into correct form for tree building
input_list <- treebuilding_preprocess(input_table, prefix, out_dir)

#### =========== RUN TREE BUILDING ========= ####

# run main CONIPHER tree building function
sample_pyclone_tree <-      treebuilding_run(sample_input_list = input_list
                                                  , ccf_buffer = ccf_buffer
                                                  , pval_cutoff = pval_cutoff
                                                  , use_boot = use_boot
                                                  , merge_clusters = merge_clusters
                                                  , correct_cpn_clusters = correct_cpn_clusters
                                                  , adjust_noisy_clusters = adjust_noisy_clusters
                                                  , adjust_noisy_clusters_prop = adjust_noisy_clusters_prop
                                                  , min_ccf = min_ccf
                                                  , min_cluster_size = min_cluster_size
                                                  , plotting = TRUE
                                                  , run.multi.trees = run.multi.trees
)

#### =========== SAVE OUTPUT ========= ####

# Save all tree building output
if(!is.na(sample_pyclone_tree$graph_pyclone[1]))
  cat('\nSaving all treebuilding output\n')
{
  ### Plotting tree
  treebuilding_plot(sample_pyclone_tree)

  ### Creating human readable format
  ### writing all trees
  treeFile <- paste0(sample_pyclone_tree$parameters$generalSave, "allTrees.txt")
  if ("alt_trees" %in% names(sample_pyclone_tree$graph_pyclone)) {
    write.table(paste0("### ", length(sample_pyclone_tree$graph_pyclone$alt_trees), " trees"), file = treeFile, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    tmp <- sapply(seq(1, length(sample_pyclone_tree$graph_pyclone$alt_trees)), function(x) {
        write.table(paste0("# tree ", x), file = treeFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")  
        write.table(sample_pyclone_tree$graph_pyclone$alt_trees[[x]], file = treeFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    })
  } else {
    write.table(paste0("### ", 1, " trees"), file = treeFile, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    write.table(paste0("# tree ", 1), file = treeFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")  
    write.table(sample_pyclone_tree$graph_pyclone$Corrected_tree, file = treeFile, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }

  ### writing consensus branches
  consensusBranchesFile <- paste0(sample_pyclone_tree$parameters$generalSave, "consensusBranches.txt")
  write.table(Reduce(rbind, strsplit(sample_pyclone_tree$graph_pyclone$consensus_branches, split = ":")), file = consensusBranchesFile, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

  ### writing consensus relationships
  consensusRelationshipsFile <- paste0(sample_pyclone_tree$parameters$generalSave, "consensusRelationships.txt")
  write.table(Reduce(rbind, strsplit(sample_pyclone_tree$graph_pyclone$consensus_relationships, split = ":")), file = consensusRelationshipsFile, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  ### writing cluster information
  clusterInfoFile <- paste0(sample_pyclone_tree$parameters$generalSave, "clusterInfo.txt")

  clusterInfoDF <- data.frame(clusterID = names(sample_pyclone_tree$graph_pyclone$edgelength), stringsAsFactors = FALSE)
  clusterInfoDF$truncal <- ifelse(clusterInfoDF$clusterID %in% sample_pyclone_tree$graph_pyclone$trunk, TRUE, FALSE)
  clusterInfoDF$treeClust <- ifelse(clusterInfoDF$clusterID %in% unique(c(sample_pyclone_tree$graph_pyclone$Corrected_tree)), TRUE, FALSE)
  clusterInfoDF$cpnRemClust <- ifelse(clusterInfoDF$clusterID %in% sample_pyclone_tree$cpn_removed_clusters, TRUE, FALSE)
  clusterInfoDF$nMuts <- as.numeric(sample_pyclone_tree$graph_pyclone$edgelength)

  clusterInfoDF <- clusterInfoDF %>% full_join(data.frame(sample_pyclone_tree$nested_pyclone$ccf_cluster_table, stringsAsFactors = FALSE) %>% mutate(clusterID = rownames(.)) %>% pivot_longer(!clusterID, names_to = "Region", values_to = "meanCCF"), by = c("clusterID"))
  clusterInfoDF <- clusterInfoDF %>% full_join(data.frame(sample_pyclone_tree$nested_pyclone$ccf_ci_lower, stringsAsFactors = FALSE) %>% mutate(clusterID = rownames(.)) %>% pivot_longer(!clusterID, names_to = "Region", values_to = "CCF_CI_low"), by = c("clusterID", "Region"))
  clusterInfoDF <- clusterInfoDF %>% full_join(data.frame(sample_pyclone_tree$nested_pyclone$ccf_ci_upper, stringsAsFactors = FALSE) %>% mutate(clusterID = rownames(.)) %>% pivot_longer(!clusterID, names_to = "Region", values_to = "CCF_CI_high"), by = c("clusterID", "Region"))
  clusterInfoDF <- clusterInfoDF %>% full_join(data.frame(sample_pyclone_tree$clonality_out$clonality_table_corrected, stringsAsFactors = FALSE) %>% mutate(clusterID = rownames(.)) %>% pivot_longer(!clusterID, names_to = "Region", values_to = "clonality"), by = c("clusterID", "Region"))
  clusterInfoDF <- clusterInfoDF %>% full_join(data.frame(sample_pyclone_tree$clone_proportion_out$clone_proportion_table, stringsAsFactors = FALSE) %>% mutate(clusterID = rownames(.)) %>% pivot_longer(!clusterID, names_to = "Region", values_to = "clone_proportions_default"), by = c("clusterID", "Region"))
  write.table(clusterInfoDF, file = clusterInfoFile, row.names = FALSE, quote = FALSE, sep = "\t")

  ### writing clone proportion information
  cloneproportionInfoFile <- paste0(sample_pyclone_tree$parameters$generalSave, "cloneProportionsMinErrorTrees.txt")

  cp_min_sce_trees <- sample_pyclone_tree$clone_proportion_out$clone_proportions_min_sce_trees
  cloneproportionInfoList <- lapply(seq(cp_min_sce_trees), function(i){
    tree_id <- names(cp_min_sce_trees)[i]
    cp_table <- data.frame(cp_min_sce_trees[[i]], stringsAsFactors = FALSE)
    cp_table$clusterID <- rownames(cp_table)
    cp_table$treeID <- tree_id
    return(cp_table)
  })
  cloneproportionInfoDF <- do.call(rbind, cloneproportionInfoList)
  write.table(cloneproportionInfoDF, file = cloneproportionInfoFile, row.names = FALSE, quote = FALSE, sep = "\t")

  ### writing output muttable - similar to input
  input_table <- input_table %>% rename(originalCLUSTER = CLUSTER)
  if (is.null(nrow(sample_pyclone_tree$merge_clusters))) {
    input_table <- input_table %>% mutate(treeCLUSTER = originalCLUSTER)
  } else {
    input_table <- input_table %>% mutate(treeCLUSTER = originalCLUSTER)
    for (i in 1:nrow(sample_pyclone_tree$merged_clusters)) {
      input_table$treeCLUSTER <- gsub(sample_pyclone_tree$merged_clusters[i, 1], sample_pyclone_tree$merged_clusters[i, 3], input_table$treeCLUSTER)
    }
  }
  write.table(input_table, file = paste0(sample_pyclone_tree$parameters$generalSave, "treeTable.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")

  ### writing alternative trees summary metrics
  altTreeInfoFile <- paste0(sample_pyclone_tree$parameters$generalSave, "alternativeTreeMetrics.txt")

  altTreeInfoDF <- data.frame(treeID = seq(sample_pyclone_tree$graph_pyclone$alt_trees), stringsAsFactors = FALSE)
  
  altTreeInfoDF$sum_condition_error <- sapply(altTreeInfoDF$treeID, function(i) sample_pyclone_tree$graph_pyclone$alt_trees_sum_condition_error[i])
  altTreeInfoDF$SCE_ranking <- match(altTreeInfoDF$sum_condition_error, sort(unique(altTreeInfoDF$sum_condition_error)))
  altTreeInfoDF$lowest_SCE <- ifelse(altTreeInfoDF$sum_condition_error == min(altTreeInfoDF$sum_condition_error), 'Lowest SCE tree', 'Alternative tree')

  altTreeInfoDF$edge_probability_score <- sapply(altTreeInfoDF$treeID, function(i) sample_pyclone_tree$graph_pyclone$alt_trees_edge_probability[i])
  altTreeInfoDF$edge_probability_ranking <- match(altTreeInfoDF$edge_probability_score, rev(sort(unique(altTreeInfoDF$edge_probability_score))))
  altTreeInfoDF$highest_edge_probability <- ifelse(altTreeInfoDF$edge_probability_score == max(altTreeInfoDF$edge_probability_score), 'Highest edge probability tree', 'Alternative tree')
  write.table(altTreeInfoDF, file = altTreeInfoFile, row.names = FALSE, quote = FALSE, sep = "\t")

}