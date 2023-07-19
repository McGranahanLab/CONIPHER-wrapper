suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(CONIPHER))

option_list <- list(
    make_option("--input_tsv_loc", type = "character", default = NULL,
                help = "File path to input mutation table in long format", metavar = "character"),
    make_option("--out_dir", type = "character", default = NULL,
                help = "Working directory where output should be saved", metavar = "character"),
    make_option("--prefix", type = "character", default = NULL,
                help = "Sample prefix", metavar = "character"),

    make_option("--ccf_buffer", type = "integer", default = 10, 
                help = "Buffer used for CCF calculations", metavar = "character"),
    make_option("--pval_cutoff", type = "double", default = 0.01, 
                help = "P-value for copy number calculation", metavar = "character"),
    
    make_option("--use_boot", type = "logical", default = TRUE,
                help = "Should bootstrapping be used", metavar = "character"),
    make_option("--merge_clusters", type = "logical", default = TRUE, 
                help = "Should similar clusters be merged if possible", metavar = "character"),  
    make_option("--correct_cpn_clusters", type = "logical", default = TRUE, 
                help = "Should clusters driven by copy number be removed", metavar = "character"),  
    make_option("--adjust_noisy_clusters", type = "logical", default = FALSE, 
                help = "Should noisy clusters be adjusted", metavar = "character"),  
    make_option("--adjust_noisy_clusters_prop", type = "double", default = 0.05, 
                help = "What is minimum proportion of mutations should be present in a region to avoid cluster adjustment", metavar = "character"),  

    make_option("--min_ccf", type = "double", default = 0.01, 
                help = "Minimum CCF to consider a mutation as present", metavar = "character"),
    make_option("--min_cluster_size", type = "integer", default = 5, 
                help = "Minimum number of mutations in a cluster for it to be included in analysis", metavar = "character"),
    make_option("--multi_trees", type = "logical", default = TRUE, 
                help = "Should alternative trees be explored", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

conipher_treebuilding(input_tsv_loc = opt$input_tsv_loc,
                      out_dir = opt$out_dir,
                      prefix = opt$prefix,
                      ccf_buffer = opt$ccf_buffer,
                      pval_cutoff = opt$pval_cutoff,
                      use_boot = opt$use_boot,
                      merge_clusters = opt$merge_clusters,
                      correct_cpn_clusters = opt$correct_cpn_clusters,
                      adjust_noisy_clusters = opt$adjust_noisy_clusters,
                      adjust_noisy_clusters_prop = opt$adjust_noisy_clusters_prop,
                      min_ccf = opt$min_ccf,
                      min_cluster_size = opt$min_cluster_size,
                      multi_trees = opt$multi_trees)
