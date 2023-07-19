suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(CONIPHER))

option_list <- list(
    make_option("--case_id", type = "character", default = NULL,
                help = "Tumour ID", metavar = "character"),
    make_option("--prefix", type = "character", default = NULL,
                help = "Sample prefix", metavar = "character"),

    make_option("--out_dir", type = "character", default = NULL,
                help = "Working directory where output should be saved", metavar = "character"),
    make_option("--input_tsv_loc", type = "character", default = NULL,
                help = "File path to input mutation table in long format", metavar = "character"),
    make_option("--input_seg_tsv_loc", type = "character", default = NULL,
                help = "File path to input segment table used for plotting", metavar = "character"),


    make_option("--subclonal_copy_correction", type = "logical", default = TRUE,
                help = "Should subclonal copy number correction be used", metavar = "character"),
    make_option("--only_truncal_subclonal_copy_correction", type="logical", default = TRUE, 
                help = "Should only truncal subclonal copy number correction be used", metavar = "character"),  

    make_option("--pyclone_yaml_loc", type = "character", default = NULL, 
                help = "Location to a template yaml file. If null package default is used", metavar = "character"),

    make_option("--min_cluster_size", type = "integer", default = 5, 
                help = "Minimum number of mutations in a cluster to be considered", metavar = "character"),
    make_option("--multiple_test_correction", type = "logical", default = TRUE, 
                help = "Should multiple testing correction be applied for the copy number correcting mutations", metavar = "character"),
    make_option("--clean_clusters", type = "logical", default = TRUE, 
                help = "Should clusters be cleaned and merged", metavar = "character"),
    make_option("--clonal_cutOff", type = "double", default = 0.9, 
                help = "Lower threshold CCF to be considered clonal", metavar = "character"),
    make_option("--propClonal_threshold", type = "double", default = 0.25, 
                help = "Proportion of cluster that needs to be considered clonal to merge", metavar = "character"),
    
    make_option("--fix_absentCCFs", type = "logical", default = TRUE, 
                help = "Should CCF of absent mutations be set to zero", metavar = "character"),
    make_option("--driver_filter", type = "character", default = "1A,1,2A", 
                help = "What filter to use for drivers", metavar = "character"),
    make_option("--burn_in", type = "integer", default = 1000, 
                help = "Burn-in for DP clustering", metavar = "character"),
    make_option("--seed", type = "integer", default = 1024, 
                help = "Seed for pyclone", metavar = "character"),
    make_option("--nProcs", type = "integer", default = 1,
                help = "Number of cores allocated to run script in parallel", metavar = "character"),



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
    make_option("--multi_trees", type = "logical", default = TRUE, 
                help = "Should alternative trees be explored", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

conipher_run(case_id = opt$case_id,
             prefix = opt$prefix,
             out_dir = opt$out_dir,
             input_tsv_loc = opt$input_tsv_loc,
             input_seg_tsv_loc = opt$input_seg_tsv_loc,
             subclonal_copy_correction = opt$subclonal_copy_correction,
             only_truncal_subclonal_copy_correction = opt$only_truncal_subclonal_copy_corection,
             pyclone_yaml_loc = opt$pyclone_yaml_loc,
             min_cluster_size = opt$min_cluster_size,
             multiple_test_correction = opt$multiple_test_correction,
             clean_clusters = opt$clean_clusters,
             clonal_cutOff = opt$clonal_cutOff,
             propClonal_threshold = opt$propClonal_threshold,
             fix_absentCCFs = opt$fix_absentCCFs,
             driver_filter = opt$driver_filter,
             burn_in = opt$burn_in,
             seed = opt$seed,
             nProcs = opt$nProcs,
             ccf_buffer = opt$ccf_buffer,
             pval_cutoff = opt$pval_cutoff,
             use_boot = opt$use_boot,
             merge_clusters = opt$merge_clusters,
             correct_cpn_clusters = opt$correct_cpn_clusters,
             adjust_noisy_clusters = opt$adjust_noisy_clusters,
             adjust_noisy_clusters_prop = opt$adjust_noisy_clusters_prop,
             min_ccf = opt$min_ccf,
             multi_trees = opt$multi_trees)
