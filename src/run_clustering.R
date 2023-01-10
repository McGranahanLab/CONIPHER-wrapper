######################################################################################################################
# This script is designed to estimate the cancer cell fraction using PyClone
######################################################################################################################

cat('\n======== RUNNING PYCLONE SIMPLE CLUSTERING ===========\n')

####################################################################################################### Load libraries
######################################################################################################################
cat('\nLoading required packages\n')

suppressPackageStartupMessages(library(compiler))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
options(expressions = 10000)

######################################################################################################## Optpars setup
######################################################################################################################
option_list = list(
  # required input:
  make_option(c("-c", "--case_id"), type="character", default=NULL, 
              help="Tumour case ID", metavar="character"),

  make_option(c("-w", "--working_dir"), type="character", default=NULL, 
              help="Working directory", metavar="character"),
  make_option(c("--script_dir"), type="character", default=NULL,
              help="Location of folder with pyclone simple clustering helper functions", metavar = "character"),

  make_option(c("-t", "--input_tsv"), type="character", default=NULL, 
              help="A mutation table in long format (mutations x tumour regions)", metavar="character"),
  
  make_option(c("--input_seg_tsv"), type="character", default=NULL, 
              help="A copy number segment .tsv used for plotting only", metavar="character"),
  
  make_option(c("--subclonal_copy_correction"), type="logical", default="TRUE", 
              help="should subclonal copy number correction be used", metavar="character"),
  
  make_option(c("--only_truncal_subclonal_copy_correction"), type="logical", default="TRUE", 
              help="should only truncal subclonal copy number correction be used", metavar="character"),

  make_option(c("--pyclone_ccf"), type="character", default="pyclone", 
              help="if pyclone original is used, should the phyloCCF be calculated, or simply the pyclone version [pyclone or phyloCCF]", metavar="character"),
  

  make_option(c("-y", "--pyclone_yaml"), type="character", default="template.config.yaml", 
              help="A template yaml file for pyclone", metavar="character"),

  make_option(c("--min_cluster_size"), type="integer", default=5, 
              help="minimum number of mutations in a cluster to be considered", metavar="character"),
  
  make_option(c("--multiple_test_correction"), type="logical", default=TRUE, 
              help="should multiple testing correction be applied for the copy number correcting mutations", metavar="character"),
  
  make_option(c("--clean_clusters"), type="logical", default=TRUE, 
              help="should clusters be cleaned and merged?", metavar="character"),
  make_option(c("--clonal_cutOff"), type="numeric", default=0.9, 
              help="lower threshold CCF to be considered clonal", metavar="character"),
  make_option(c("--propClonal_threshold"), type="numeric", default=0.25, 
              help="proportion of cluster that needs to be considered clonal to merge", metavar="character"),


  make_option(c("--fix_absentCCFs"), type="logical", default=TRUE, 
              help="should CCF of absent mutations be set to zero?", metavar="character"),

  
  make_option(c("--driver_filter"), type="character", default="1A,1,2A", 
              help="what filter to use for drivers", metavar="character"),
  make_option(c("--burn_in"), type="integer", default=1000, 
              help="burn-in for DP clustering", metavar="character"),
  make_option(c("--nProcs"), type = "integer", default = 1,
              help="Number of cores allocated to run script in parallel", metavar="character")
)

opt_parser       <- OptionParser(option_list=option_list)
opt              <- parse_args(opt_parser)

print(opt)


cat(opt$input_tsv)
# DIRECTORIES ####
cat('\nSetting up working Directories')
if (!file.exists(opt$working_dir)) {
  dir.create(opt$working_dir, recursive = TRUE)
}
setwd(opt$working_dir)
# make sure working directory ends in a '/'
if (substr(opt$working_dir, nchar(opt$working_dir), nchar(opt$working_dir)) != '/') opt$working_dir <- paste0(opt$working_dir, '/')


# PARAMETER SETTING ####
cat('\nSetting parameters')
# set parameters used for script (these are currently hardcoded)

template.config.yaml <- paste0(opt$script_dir, "/", opt$pyclone_yaml)
driver.cat           <- unlist(strsplit(opt$driver_filter,split=","))
min.cluster.size     <- opt$min_cluster_size
burn_in              <- opt$burn_in
new.dir              <- opt$working_dir
patient              <- opt$case_id
sample               <- patient
skip.VAFQC           <- TRUE
gender               <- 'male' # at the moment gender isn't really used, so hardcoding as male
min.ccf.present      <- 0.1
PyClone              <- "PyClone"
run.pyclone          <- TRUE
nProcs               <- opt$nProcs
scriptDir            <- opt$script_dir


source(paste0(scriptDir, "TRACERxHelperFunctions.R"))
require(KernSmooth)



###### Run script #####################


cat('\nPyclone simple clustering analysis of the following tumour case:\n')
print(patient)
cat('\n')


if( !file.exists(new.dir))
{
  if( !dir.create(new.dir, recursive = TRUE) )
  {
    stop("Unable to create root directory.\n")
  }
}


input_tsv_loc <- opt$input_tsv

cat(opt$working_dir)


if( !file.exists(input_tsv_loc))
{
  cat(input_tsv_loc)
  stop("Unable to find mut.table.loc.\n")
  
}

input_tsv <- read.delim(input_tsv_loc, sep = "\t", stringsAsFactors = FALSE, header = TRUE, fill = TRUE, quote = "")
input_tsv[, "key"] <- paste(paste0("chr", input_tsv[, "CHR"]), input_tsv[, "POS"], input_tsv[, "REF"], input_tsv[, "ALT"], sep = ":")

if (sum(grepl("Use.For.Plots", colnames(input_tsv))) == 2) {
  ### create input needed for pyclone
  mut.table <- data.frame(key = input_tsv[, "key"],
                          chr = input_tsv[, "CHR"],
                          start = input_tsv[, "POS"],
                          stop = input_tsv[, "POS"],
                          ref = input_tsv[, "REF"],
                          var = input_tsv[, "ALT"],
                          is_SNV = TRUE,
                          Use.For.Plots = input_tsv[, "Use.For.Plots"],
                          Use.For.Plots.Indel = input_tsv[, "Use.For.Plots.Indel"],
                          stringsAsFactors = FALSE)

  mut.table <- mut.table %>% 
    full_join(input_tsv %>% 
      select(key, SAMPLE, REF_COUNT, VAR_COUNT, DEPTH) %>% 
      dplyr::rename(cov = DEPTH, ref_count = REF_COUNT, var_count = VAR_COUNT) %>% 
      dplyr::mutate(VAF = var_count / cov * 100) %>%
        tidyr::pivot_wider(names_from = SAMPLE, values_from = c(cov, ref_count, var_count, VAF), names_glue = "{SAMPLE}.{.value}"), by = "key") %>%
    rowwise() %>%
    mutate(is_blacklist = FALSE,
      max.var_count = max(dplyr::c_across(ends_with(".var_count"))),
      max.VAF = max(dplyr::c_across(ends_with(".VAF"))))

  mut.table <- mut.table %>%
    mutate(Gene.refGene = "") %>%
    mutate(driverCategory = NA)

  mut.table <- data.frame(mut.table, stringsAsFactors = FALSE)
} else {
  ### create input needed for pyclone
  mut.table <- data.frame(key = input_tsv[, "key"],
                          chr = input_tsv[, "CHR"],
                          start = input_tsv[, "POS"],
                          stop = input_tsv[, "POS"],
                          ref = input_tsv[, "REF"],
                          var = input_tsv[, "ALT"],
                          is_SNV = TRUE,
                          stringsAsFactors = FALSE)

  mut.table <- mut.table %>% 
    full_join(input_tsv %>% 
      select(key, SAMPLE, REF_COUNT, VAR_COUNT, DEPTH) %>% 
      dplyr::rename(cov = DEPTH, ref_count = REF_COUNT, var_count = VAR_COUNT) %>% 
      dplyr::mutate(VAF = var_count / cov * 100) %>%
        tidyr::pivot_wider(names_from = SAMPLE, values_from = c(cov, ref_count, var_count, VAF), names_glue = "{SAMPLE}.{.value}"), by = "key") %>%
    rowwise() %>%
    mutate(is_blacklist = FALSE,
      max.var_count = max(dplyr::c_across(ends_with(".var_count"))),
      max.VAF = max(dplyr::c_across(ends_with(".VAF"))))

  mut.table <- mut.table %>%
    mutate(Gene.refGene = "") %>%
    mutate(Use.For.Plots = TRUE, Use.For.Plots.Indel = FALSE) %>%
    mutate(driverCategory = NA)

  mut.table <- data.frame(mut.table, stringsAsFactors = FALSE)
}

print('test6')

if(TRUE%in%grepl('CN_A',colnames(input_tsv))) # KG:  delete - this shouldnt be needed if we have provided correct input columns.
{
  colnames(input_tsv) <- gsub("CN_A","COPY_NUMBER_A",colnames(input_tsv))
  colnames(input_tsv) <- gsub("CN_B","COPY_NUMBER_B",colnames(input_tsv))
  colnames(input_tsv) <- gsub("PURITY","ACF",colnames(input_tsv))
  
}

seg.mat.copy <- data.frame(SampleID = input_tsv[,"SAMPLE"],
                           chr = input_tsv[, "CHR"],
                           startpos = input_tsv[, "POS"],
                           endpos = input_tsv[, "POS"],
                           n.het = 1,
                           cnTotal = round(input_tsv[, "COPY_NUMBER_A"] + input_tsv[, "COPY_NUMBER_B"]),
                           nMajor = round(input_tsv[, "COPY_NUMBER_A"]),
                           nMinor = round(input_tsv[, "COPY_NUMBER_B"]),
                           Ploidy = input_tsv[, "PLOIDY"],
                           ACF = input_tsv[, "ACF"],
                           COPY_NUMBER_A = input_tsv[, "COPY_NUMBER_A"],
                           COPY_NUMBER_B = input_tsv[, "COPY_NUMBER_B"], 
                           stringsAsFactors = FALSE)



regions <- unique(seg.mat.copy$SampleID)
regions <- gsub("-","\\.",regions)
print(regions)

        mut.table      <- mut.table[mut.table$chr%in%1:22,,drop=FALSE]
        seg.mat.copy   <- seg.mat.copy[seg.mat.copy$chr%in%1:22,,drop=FALSE]
pat<-patient


  regions.to.use  <- unique(seg.mat.copy$Sample)

mut.table$mutation_id  <- paste(sample,mut.table$chr,mut.table$start,mut.table$ref,sep=":")
mut.table              <- mut.table[order(mut.table$max.VAF,decreasing=TRUE),]
mut.table              <- mut.table[!duplicated(mut.table$mutation_id),,drop=FALSE]
rownames(mut.table)    <- mut.table$mutation_id




max.vaf        <- c()
max.var.count  <- c()

for (region in regions.to.use)
{
  
  if(paste(region, ".VAF",sep="")%in%colnames(mut.table))
  {
    max.vaf             <- cbind(max.vaf,mut.table[,paste(region, ".VAF",sep="")])
    max.var.count       <- cbind(max.var.count,mut.table[,paste(region, ".var_count",sep="")])
    
  }
  
}

mut.table$max.VAF        <- apply(max.vaf,1,max)
mut.table$max.var_count  <- apply(max.var.count,1,max)

mut.table                <- mut.table[!is.na(mut.table$max.VAF),,drop=FALSE]

mut.table <- mut.table[mut.table$Use.For.Plots | mut.table$Use.For.Plots.Indel, ]
mut.table <- mut.table[!((mut.table$Use.For.Plots & mut.table$max.var_count < 10) | is.na(mut.table$max.var_count)), ]

seg.mat.phylo      <- create.subclonal.copy.number(seg.mat.copy = seg.mat.copy,min.subclonal = 0.01)

if(opt$subclonal_copy_correction%in%"FALSE")
{
  cat('\nRunning without subclonal copy number mode')
  seg.mat.phylo$COPY_NUMBER_A <- seg.mat.phylo$nMajor
  seg.mat.phylo$COPY_NUMBER_B <- seg.mat.phylo$nMinor
  seg.mat.phylo$fracA <- 1
  seg.mat.phylo$fracB <- 0
  seg.mat.phylo$fracC <- NA
  seg.mat.phylo$fracD <- NA
  seg.mat.phylo$nMaj_A <- seg.mat.phylo$nMajor
  seg.mat.phylo$nMin_A <- seg.mat.phylo$nMinor
  seg.mat.phylo$nMaj_B <- seg.mat.phylo$nMajor
  seg.mat.phylo$nMin_B <- seg.mat.phylo$nMinor
  seg.mat.phylo$nMaj_C <- NA
  seg.mat.phylo$nMin_C <- NA
  seg.mat.phylo$nMaj_D <- NA
  seg.mat.phylo$nMin_D <- NA
  
  seg.mat.copy$COPY_NUMBER_A <- seg.mat.phylo$nMajor
  seg.mat.copy$COPY_NUMBER_B <- seg.mat.phylo$nMinor
  
  
}


if(nrow(mut.table)==0)
{
  stop('No mutations passed filtering, stopping PyClone phylo clustering')
}

patient.list          <- list()
phylo.region.list     <- list()
cellularity           <- rep(NA,length(regions.to.use))


# Let's plot the copy number. 
patient.list <- list()
phylo.region.list <- list()
cellularity           <- rep(NA,length(regions.to.use))
names(cellularity)    <- regions.to.use

# determine the indelCorrectionFactor
if (length(regions.to.use)>1){
    indelCorrectionFactor <- determineIndelCorrectionFactor(mut.table = mut.table,regions.to.use = regions.to.use,seg.mat.phylo = seg.mat.phylo)
    indelMuts <- rownames(mut.table[mut.table$Use.For.Plots.Indel%in%TRUE,,drop=FALSE])
}

pdf(paste(new.dir, sample, ".subclonal.mut.cpn.pdf",sep=""))
lyout <- c()
for (i in seq(1,length(regions.to.use)*2,by=2))
{
  lyout <- rbind(lyout,rbind(c(rep(i,9),i+1)))
}

layout(lyout)

if(nrow(lyout)==1)
{
  par(mar=c(9,4,9,4))
}

for (region in regions.to.use)
{
  
  region.mut.table <- mut.table
  
  region.seg.copy  <- seg.mat.copy[seg.mat.copy$SampleID%in%region,,drop=FALSE]
  region.seg.phylo <- seg.mat.phylo[seg.mat.phylo$SampleID%in%region,,drop=FALSE]
  pyclone.table    <- data.frame(t(sapply(1:nrow(region.mut.table),identify.subclonal.mut.copy.number.ascat,region.mut.table,region.seg.phylo,region,sample))
                                 ,stringsAsFactors=FALSE)
  pyclone.table    <- pyclone.table[!is.na(pyclone.table$minor_cn),]
  pyclone.table    <- pyclone.table[!is.na(pyclone.table$ref_counts),]
  pyclone.table    <- pyclone.table[!duplicated(pyclone.table$mutation_id),]
  
       sample.purity <- region.seg.copy$ACF[1]

  pyclone.table              <- pyclone.table[(as.numeric(pyclone.table$ref_counts)+as.numeric(pyclone.table$var_counts))>=1,,drop=FALSE]
  region.earlyLate           <- earlyORlateGender(region=region,complete.mutation.table=pyclone.table,purity=sample.purity,gender=gender)
  if(opt$multiple_test_correction%in%FALSE)
  {
    region.phyloCCF            <- calculate_phylo_ccf(region=region
                                                            ,complete.mutation.table=pyclone.table
                                                            ,purity=sample.purity
                                                            ,order.by.pos = TRUE
                                                            ,gender=gender
    )
  }

  if(opt$multiple_test_correction%in%TRUE)
  {
    region.phyloCCF            <- calculate_phylo_ccf_withBH(region=region
                                                            ,complete.mutation.table=pyclone.table
                                                            ,purity=sample.purity
                                                            ,order.by.pos = TRUE
                                                            ,gender=gender
    )
  }
  
  if (length(regions.to.use)>1){
    tmp <- intersect(rownames(region.phyloCCF),indelMuts)
     # let's look at indels specifically
     region.phyloCCF[tmp,'phyloCCF'] <- region.phyloCCF[tmp,'phyloCCF']*as.numeric(indelCorrectionFactor[region])
     region.phyloCCF[tmp,'mutCopyNum'] <- region.phyloCCF[tmp,'mutCopyNum']*as.numeric(indelCorrectionFactor[region])
  }
  
  phylo.region.list[[region]]                                              <- region.phyloCCF
  
  cellularity[region]                                                      <- sample.purity
  
  # add extra column to 
  sub.mat.copy               <- region.seg.copy
  sub.mat.copy$Copy.Number   <- sub.mat.copy$cnTotal
  sub.mat.copy$min.allele    <- apply(cbind(sub.mat.copy$nMinor,sub.mat.copy$nMajor),1,min)
  
  
  colnames(sub.mat.copy)[2]  <- 'Chromosome'
  colnames(sub.mat.copy)[3]  <- 'StartPosition'
  colnames(sub.mat.copy)[4]  <- 'EndPosition'
  colnames(sub.mat.copy)[5]  <- 'nr.probes'
  
  
  par(mar=c(0.5,4,0.5,0.2))
  if(nrow(lyout)==1)
  {
    par(mar=c(9,4,9,0.2))
  }
  
  par(lend=1)
  plot.EarlyOrLate.raw(seg.mat.patient=sub.mat.copy
                       ,TCGA.earlyLate=region.phyloCCF
                       ,TCGA.purity=sample.purity
                       ,sub.clonal=1
  )
  mtext(unlist(strsplit(region,split="\\."))[2], side=2,cex=1,line=2,las=2)
  
  ds <- density(ifelse(as.numeric(region.phyloCCF$mutCopyNum)>5,5,as.numeric(region.phyloCCF$mutCopyNum)))
  ds1 <- ds
  ds1$x <- ds$y
  ds1$y <- ds$x
  par(mar=c(0.5,0,0.5,4))
  if(nrow(lyout)==1)
  {
    par(mar=c(9,0,9,4))
  }
  A <- hist(ifelse(as.numeric(region.phyloCCF$mutCopyNum)>5,5,as.numeric(region.phyloCCF$mutCopyNum)),breaks=seq(-0.25,6,by=0.1),plot=FALSE)
  plot(NULL, type = "n"
       , xlim = c(0, max(A$density))
       , ylim=c(-0.25,6)
       ,bty='n'
       ,xaxs='i'
       ,xaxt='n'
       ,yaxt='n'
       ,yaxs='i'
       ,xlab=""
       ,main=""
       ,ylab=""
  )
  rect(0, A$breaks[1:(length(A$breaks) - 1)], A$density, A$breaks[2:length(A$breaks)]
       ,border=TRUE,col="#CC6666")
  
  lines(ds1)
  
}
dev.off()


# perform additional copy number correction where needed
muts_to_consider <- unlist(phylo.region.list[[1]]$mutation_id)
for (mut in muts_to_consider)
{
    mut_table <- c()
    for (region in names(phylo.region.list))
    {
      mut_table <- rbind(mut_table,phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%mut,,drop=FALSE])
    }
    
    if(max(mut_table$phyloCCF)>1.5)
    {
      #all mutations are greater than 1.5 --> re-centre on 1 or       #one or more (but not all) mutation is greater than 1.5 (but all are clonal) --> re-centre on 1
      if(length(which(mut_table$phyloCCF>1.5))==nrow(mut_table)|length(which(mut_table$phyloCCF.0.05>=1))==nrow(mut_table))
      {
        small_mut_table <- mut_table[mut_table$phyloCCF.0.05>1,,drop=FALSE]
        if(nrow(small_mut_table)==0)
        {next;}
        for (i in 1:nrow(small_mut_table))
        {
          
          small_row <- small_mut_table[i,,drop=FALSE]
          region <- unlist(small_row$region)
          region.copy <- seg.mat.phylo[seg.mat.phylo$SampleID%in%region,]
          
          phyloCCF      <- small_row$absolute.ccf
          phyloCCF_0.05 <- min(small_row$absolute.ccf.0.05,small_row$absolute.ccf.0.05-abs(small_row$phyloCCF-small_row$phyloCCF.0.05))
          phyloCCF_0.95 <- max(small_row$absolute.ccf.0.95,small_row$absolute.ccf.0.95+abs(small_row$phyloCCF-small_row$phyloCCF.0.05))
          
          phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%small_row$mutation_id,]$phyloCCF.0.05 <- phyloCCF_0.05
          phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%small_row$mutation_id,]$phyloCCF.0.95 <- phyloCCF_0.95
          phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%small_row$mutation_id,]$phyloCCF <- phyloCCF 
        } 
      }
      
      #at least one mutation is subclonal --> properly adjust no.chrs.bearing.mut
      if(length(which(mut_table$phyloCCF.0.95<1))>=1)
      {
        small_mut_table <- mut_table[mut_table$phyloCCF.0.05>1,,drop=FALSE]
        if(nrow(small_mut_table)==0)
        {next;}
        for (i in 1:nrow(small_mut_table))
        {
          
          small_row <- small_mut_table[i,,drop=FALSE]
          region <- unlist(small_row$region)
          region.copy <- seg.mat.phylo[seg.mat.phylo$SampleID%in%region,]
          
          phyloCCF      <- small_row$phyloCCF/small_row$no.chrs.bearing.mut
          phyloCCF_0.05 <- small_row$phyloCCF.0.05/small_row$no.chrs.bearing.mut
          phyloCCF_0.95 <- small_row$phyloCCF.0.95/small_row$no.chrs.bearing.mut
          phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%small_row$mutation_id,]$phyloCCF.0.05 <- phyloCCF_0.05
          phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%small_row$mutation_id,]$phyloCCF.0.95 <- phyloCCF_0.95
          phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%small_row$mutation_id,]$phyloCCF <- phyloCCF
          
        }    
      } 
    }
    
    if(!TRUE%in%c(mut_table$phyloCCF>mut_table$mutCopyNum))
    {
     next;
    }
    if(TRUE%in%c(mut_table$phyloCCF>mut_table$mutCopyNum))
    {
      
      # are all the mutations now truncal?
      if(length(which(mut_table$phyloCCF.0.95<1))==0)
      {
        next;
      }
      if(opt$only_truncal_subclonal_copy_correction%in%TRUE)
      {
        if(length(which(mut_table$phyloCCF.0.95<1))>=1)
        {
        
          
          # so we've performed copy number correction, but that didn't make the mutation clonal, so let's revert back to non-copy number corrected
          muts_to_revert <- mut_table[mut_table$no.chrs.bearing.mut<1,,drop=FALSE]
          
          for (region in unlist(muts_to_revert$region))
          {
            mut_to_revert <- muts_to_revert[unlist(muts_to_revert$region)%in%region,,drop=FALSE]
            region.copy <- seg.mat.phylo[seg.mat.phylo$SampleID%in%region,]
            expVAF   <- min(1 - 1e-6, c((region.copy$ACF[1]*1) / (2*(1-region.copy$ACF[1]) + region.copy$ACF[1]*(as.numeric(mut_to_revert$major_raw)+as.numeric(mut_to_revert$minor_raw)))))
            VAF_ci   <- prop.test(x = as.numeric(mut_to_revert$var_counts),n = as.numeric(mut_to_revert$ref_counts)+as.numeric(mut_to_revert$var_counts),p = expVAF)
            phyloCCF      <- (VAF_ci$estimate *1/region.copy$ACF[1])*((region.copy$ACF[1]*(as.numeric(mut_to_revert$major_raw)+as.numeric(mut_to_revert$minor_raw)))+2*(1-region.copy$ACF[1]))
            phyloCCF_0.05 <- (VAF_ci$conf.int[1] *1/region.copy$ACF[1])*((region.copy$ACF[1]*(as.numeric(mut_to_revert$major_raw)+as.numeric(mut_to_revert$minor_raw)))+2*(1-region.copy$ACF[1]))
            phyloCCF_0.95 <- (VAF_ci$conf.int[2] *1/region.copy$ACF[1])*((region.copy$ACF[1]*(as.numeric(mut_to_revert$major_raw)+as.numeric(mut_to_revert$minor_raw)))+2*(1-region.copy$ACF[1]))
          
            phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%mut_to_revert$mutation_id,]$no.chrs.bearing.mut <- 1
            phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%mut_to_revert$mutation_id,]$phyloCCF.0.05 <- phyloCCF_0.05
            phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%mut_to_revert$mutation_id,]$phyloCCF.0.95 <- phyloCCF_0.95
            phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%mut_to_revert$mutation_id,]$phyloCCF <- phyloCCF
            phylo.region.list[[region]][phylo.region.list[[region]]$mutation_id%in%mut_to_revert$mutation_id,]$expected.VAF <- expVAF
            
          }
          
          
          
        }
      }
      if(opt$only_truncal_subclonal_copy_correction%in%FALSE)
      {
        next;
      }
    }
}


mut.heatmap                <- matrix(NA,nrow=nrow(mut.table),ncol=length(regions.to.use))
rownames(mut.heatmap)      <- mut.table$mutation_id
colnames(mut.heatmap)      <- regions.to.use


for (region in regions.to.use)
{
  present.absent        <- unlist(phylo.region.list[[region]]$absolute.ccf)
  names(present.absent) <- phylo.region.list[[region]]$mutation_id
  present.absent        <- ifelse(present.absent>min.ccf.present,1,0)
  present.absent        <- present.absent[names(present.absent)%in%rownames(mut.heatmap)]
  mut.heatmap[names(present.absent),region] <- present.absent
  
}

mut.heatmap                <- mut.heatmap[!is.na(rowSums(mut.heatmap)),,drop=FALSE]


#, let's see if we can produce a better heatmap, with the clusters
clusters             <- apply(mut.heatmap,1,PasteVector,sep=":")

if (length(unique(apply(mut.heatmap,1,PasteVector,sep=":")))==1)
{
  warning("Only one cluster with size greater than min.cluster.size\n")
}

if(ncol(mut.heatmap)==2)
{
  #let's reduce the complexity
  AllIds <- apply(mut.heatmap,1,PasteVector,sep=":")

  if (length(unique(apply(mut.heatmap,1,PasteVector,sep=":")))==1) {
    SmallMutHeatmap <- mut.heatmap[!duplicated(AllIds),, drop = FALSE]
  } else {
    SmallMutHeatmap <- mut.heatmap[!duplicated(AllIds),]
    SmallHt <- heatmap.2(apply(SmallMutHeatmap,2,as.numeric),col=c('gray','blue'))
  }  

  RowInd <- c()
  k <- 1
  if (length(unique(apply(mut.heatmap,1,PasteVector,sep=":")))==1) {
    RowInd <- c(RowInd, which(rownames(mut.heatmap) %in% names(AllIds)))
    k <- k + 1
  }
  while(k <= length(unique(AllIds)))
  {
    RowInd <- c(RowInd,which(rownames(mut.heatmap)%in%names(AllIds[AllIds%in%names(table(AllIds)[c(AllIds[!duplicated(AllIds)][SmallHt$rowInd])][k])])))
    k <- k+1  
  }
   
  ccf.heatmap                <- matrix(NA,nrow=nrow(mut.table),ncol=length(regions.to.use))
  rownames(ccf.heatmap)      <- mut.table$mutation_id
  colnames(ccf.heatmap)      <- regions.to.use
  
  
  for (region in regions.to.use)
  {
    present.absent        <- unlist(phylo.region.list[[region]]$absolute.ccf)
    names(present.absent) <- phylo.region.list[[region]]$mutation_id
    present.absent        <- present.absent[names(present.absent)%in%rownames(mut.heatmap)]
    ccf.heatmap[names(present.absent),region] <- present.absent
    
  }
  
  ccf.heatmap                <- ccf.heatmap[!is.na(rowSums(apply(ccf.heatmap,2,as.numeric))),,drop=FALSE]
  
  clusters                   <- apply(mut.heatmap,1,PasteVector,sep=":")
  min.clusters               <- unique(apply(mut.heatmap,1,PasteVector,sep=":"))
  names(min.clusters)        <- c(1:length(min.clusters))
  most.likely.cluster        <- match(clusters,min.clusters)
  names(most.likely.cluster) <- rownames(mut.heatmap)
  
  
  
  heatmap.rows          <- rownames(mut.heatmap[rev(RowInd),])
  heatmap.rows          <- heatmap.rows[heatmap.rows%in%names(most.likely.cluster)]
  cluster.order         <- most.likely.cluster[heatmap.rows][!duplicated(most.likely.cluster[heatmap.rows])]
  heatmap.ccf           <- c()
  for (cluster in cluster.order)
  {
    ccf.ht       <- ccf.heatmap[names(most.likely.cluster[most.likely.cluster==cluster]),,drop=FALSE]
    if(nrow(ccf.ht)==1)
    {
      heatmap.ccf  <- rbind(heatmap.ccf,ccf.ht)
    }
    
    if(nrow(ccf.ht)>1)
    {
      heatmap.ccf  <- rbind(heatmap.ccf,ccf.ht[order(rowSums(apply(ccf.ht,2,as.numeric)),decreasing=TRUE),])  
    }
    
    
  }
  
}

if(ncol(mut.heatmap)>2)
{
  
  ht <- heatmap(apply(mut.heatmap,2,as.numeric),col=c('gray','blue'))
  
  
  ccf.heatmap                <- matrix(NA,nrow=nrow(mut.table),ncol=length(regions.to.use))
  rownames(ccf.heatmap)      <- mut.table$mutation_id
  colnames(ccf.heatmap)      <- regions.to.use
  
  
  for (region in regions.to.use)
  {
    present.absent        <- unlist(phylo.region.list[[region]]$absolute.ccf)
    names(present.absent) <- phylo.region.list[[region]]$mutation_id
    present.absent        <- present.absent[names(present.absent)%in%rownames(mut.heatmap)]
    ccf.heatmap[names(present.absent),region] <- present.absent
    
  }
  
  ccf.heatmap                <- ccf.heatmap[!is.na(rowSums(apply(ccf.heatmap,2,as.numeric))),,drop=FALSE]
  
  
  
  clusters                   <- apply(mut.heatmap,1,PasteVector,sep=":")
  min.clusters               <- unique(apply(mut.heatmap,1,PasteVector,sep=":"))
  names(min.clusters)        <- c(1:length(min.clusters))
  most.likely.cluster        <- match(clusters,min.clusters)
  names(most.likely.cluster) <- rownames(mut.heatmap)
  
  
  
  heatmap.rows          <- rownames(mut.heatmap[rev(ht$rowInd),])
  heatmap.rows          <- heatmap.rows[heatmap.rows%in%names(most.likely.cluster)]
  cluster.order         <- most.likely.cluster[heatmap.rows][!duplicated(most.likely.cluster[heatmap.rows])]
  heatmap.ccf           <- c()
  for (cluster in cluster.order)
  {
    ccf.ht       <- ccf.heatmap[names(most.likely.cluster[most.likely.cluster==cluster]),,drop=FALSE]
    if(nrow(ccf.ht)==1)
    {
      heatmap.ccf  <- rbind(heatmap.ccf,ccf.ht)
    }
    
    if(nrow(ccf.ht)>1)
    {
      heatmap.ccf  <- rbind(heatmap.ccf,ccf.ht[order(rowSums(apply(ccf.ht,2,as.numeric)),decreasing=TRUE),])  
    }
    
    
  }
  
}



library(parallel)
library(future)
no_cores <- nProcs
print(paste0("Number of cores that are available: ", no_cores))

source(paste0(scriptDir, "functionsForSimpleClustering.v13.R"))

simpleClusterList <- findSimpleClusters(phylo.region.list, mut.table)
### always run pyclone
tmp <- mclapply(simpleClusterList, function(x) {
  if (length(x$MutationsWithCluster) < 5) {
    CreateOutputNoPyCloneRun(clusterName = x$clusterID, patientID = sample, SmallClusters = simpleClusterList, patientDirToUse = new.dir)
  } else {
    RunPyCloneWithSimpleClusters(clusterName = x$clusterID, patientID = sample, SmallClusters = simpleClusterList, patientDirToUse = new.dir, yamlConfigLoc = template.config.yaml, run.pyclone = TRUE, pyclone.module = "PyClone/0.12.3-foss-2016b-Python-2.7.12-tkinter")
  }
}, mc.cores = no_cores)
# }
rm(list = c("no_cores", "tmp"))


### changing CCFs to 0 if cluster is absent
if (opt$fix_absentCCFs) {
  clusterPresence.df <- Reduce(rbind, lapply(names(simpleClusterList), function(x) data.frame(clusterID = x, regions = simpleClusterList[[x]]$RegionsInCluster, stringsAsFactors = FALSE)))
  for (region in names(phylo.region.list)) {
    tmp.clusterPresence <- clusterPresence.df %>% filter(regions %in% region) %>% pull(clusterID)
    tmp.clusterMutations <- as.character(unlist(lapply(simpleClusterList[as.character(tmp.clusterPresence)], function(x) x$MutationsWithCluster)))
    phylo.region.list[[region]][!as.character(unlist(phylo.region.list[[region]]$mutation_id)) %in% tmp.clusterMutations, "phyloCCF"] <- 0
    phylo.region.list[[region]][!as.character(unlist(phylo.region.list[[region]]$mutation_id)) %in% tmp.clusterMutations, "phyloCCF.0.05"] <- 0
    phylo.region.list[[region]][!as.character(unlist(phylo.region.list[[region]]$mutation_id)) %in% tmp.clusterMutations, "phyloCCF.0.95"] <- 0
  }
}
ITH1clust <- names(which(sapply(simpleClusterList, function(x) length(x$RegionsInCluster)) == length(phylo.region.list)))
ITH1muts  <- simpleClusterList[[as.character(ITH1clust)]]$MutationsWithCluster


allClusters <- paste0(list.files(new.dir, pattern = paste0(sample, "_cluster"), full.names = TRUE), "/", sample, ".results.tsv")
pyclone.results.list <- lapply(names(simpleClusterList), function(clusterID) {
  cluster.results.file <- grep(paste0("cluster", clusterID, "/"), allClusters, value = TRUE)
  pyclone.results <- read.table(cluster.results.file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  colnames(pyclone.results)  <- gsub('mutation_id', 'X', colnames(pyclone.results))
  pyclone.results$cluster_id <- as.numeric(clusterID) * 100 + as.numeric(pyclone.results$cluster_id)
  if (length(grep("^X$", colnames(pyclone.results))) == 0) {
    pyclone.results$X <- rownames(pyclone.results)
  }
  return(pyclone.results[, c("X", "cluster_id")])
})
pyclone.results <- Reduce(rbind, pyclone.results.list)
sample.results <- paste(new.dir,"/",sample,'.results.tsv',sep="")
write.table(pyclone.results, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, file = sample.results)

pyclone.results <- read.table(sample.results
                              ,sep="\t"
                              ,header=TRUE
                              ,stringsAsFactors=FALSE
)

# let's make sure the same mutations are being used. 
colnames(pyclone.results)  <- gsub('mutation_id', 'X', colnames(pyclone.results))
rownames(pyclone.results)  <- pyclone.results$X
pyclone.results            <- pyclone.results[rownames(pyclone.results)%in%mut.table$mutation_id,,drop=FALSE]

most.likely.cluster        <- pyclone.results$cluster_id
names(most.likely.cluster) <- pyclone.results$X
muts.to.remove <- c()
# let's have a look at whether any of the clusters are explained by copy number events


mut.pvals            <- c()
cluster.prop.aber    <- c()
require(coin)

for (cluster in unique(most.likely.cluster))
{
  
  tmp               <- copy.driven.clusterNEW(cluster
                                            ,seg.mat.copy=seg.mat.copy
                                            ,most.likely.cluster=most.likely.cluster
                                            ,region.earlyLate.list=phylo.region.list
                                            ,min.prop.cens=0.1
                                            ,loss.thresh = 0.25
                                            ,diff.thresh =0.55
  )
  
  if(TRUE%in%is.na(tmp[,1]))
  {
    stop;
  }
  cluster.prop.aber <- c(cluster.prop.aber,length(which(tmp[,1]<0.05&tmp[,2]<=0.5))/nrow(tmp))
  mut.pvals <- rbind(mut.pvals,tmp)
  
  
}

names(cluster.prop.aber)  <- unique(most.likely.cluster)
# Let's also check what the minimum copy number is for each mutation that is classified as lost

out <- c()
for (mutation_id in rownames(mut.pvals))
{
  out <- rbind(out,getMinCPN(mutation_id,phylo.region.list))
}

mut.pvals                 <- cbind(mut.pvals,out[,2])
sig.pvals                 <- mut.pvals[which(as.numeric(mut.pvals[,1])<0.05&as.numeric(mut.pvals[,2])<=0.75&as.numeric(mut.pvals[,3])==0),1]
sig.table                 <- mut.pvals[which(as.numeric(mut.pvals[,1])<0.05&as.numeric(mut.pvals[,2])<=0.75&as.numeric(mut.pvals[,3])==0),,drop=FALSE]

print(length(sig.pvals))
# let's see what happens if we set any cluster with over X% sig pvals as the entire cluster missing
clusters.to.remove <- names(which(as.numeric(cluster.prop.aber)>=0.85))
muts.to.remove     <- names(most.likely.cluster)[which(most.likely.cluster%in%clusters.to.remove)]
muts.to.remove     <- unique(c(names(sig.pvals),muts.to.remove))
muts.to.remove     <- muts.to.remove[!is.na(muts.to.remove)]
names(muts.to.remove) <- muts.to.remove

# let's be clever about this and give the clusters we're removing the highest names
cluster.size.remove          <- names(most.likely.cluster[most.likely.cluster%in%names(which(table(most.likely.cluster)<min.cluster.size))])
small.clusters               <- most.likely.cluster[cluster.size.remove]

# are any of the small clusters private
for(small.cluster in unique(small.clusters))
{
  
  # create a presence absence heatmap
  binary.heatmap  <- pyclone.results[,2:(ncol(pyclone.results)-1),drop=FALSE]
  binary.heatmap  <-ifelse(binary.heatmap<min.ccf.present,0,1)
  rownames(binary.heatmap) <- pyclone.results$X
  row.names       <- pyclone.results[pyclone.results$cluster_id%in%small.cluster,'X']
  
  if(median(rowSums(binary.heatmap[row.names,,drop=FALSE]))==1)
  {
    small.clusters <- small.clusters[!small.clusters%in%small.cluster]
  }
  
  
}

cluster.size.remove          <- cluster.size.remove[cluster.size.remove%in%names(small.clusters)]


old.most.likely.cluster      <- most.likely.cluster

most.likely.cluster <- most.likely.cluster[!names(most.likely.cluster)%in%unique(c(muts.to.remove,cluster.size.remove))]

if(length(most.likely.cluster)==0)
{
  stop("ClusterSize\nYou don't have any clean clusters greater than min.cluster.size")
}

# Let's rename the most likely cluster
tmp <- most.likely.cluster
most.likely.cluster <- as.numeric(factor(most.likely.cluster))
names(most.likely.cluster) <- names(tmp)
if(length(most.likely.cluster)>1)
{
  tmp                 <- as.character(1:length(unique(most.likely.cluster)))
  names(tmp)          <- names(table(most.likely.cluster)[order(table(most.likely.cluster),decreasing = TRUE)])
  tmp2                <- most.likely.cluster
  most.likely.cluster <- tmp[as.character(match(most.likely.cluster,1:length(most.likely.cluster)))]
  names(most.likely.cluster) <- names(tmp2)
  solid.cluster.end <- max(as.numeric(most.likely.cluster))
  new.cluster.start <- max(as.numeric(most.likely.cluster))+1
  removed.muts      <- old.most.likely.cluster[which(!names(old.most.likely.cluster)%in%names(most.likely.cluster))]
  removed.clusters  <- c()    
  for (cluster in names(table(removed.muts)[order(table(removed.muts),decreasing=TRUE)]))
  {
    new.cluster <- removed.muts[removed.muts%in%cluster]
    new.cluster[new.cluster%in%cluster] <- new.cluster.start
    removed.clusters <- c(removed.clusters,new.cluster)
    new.cluster.start <- new.cluster.start +1
  }
  
  # let's add back the ones to remove, but make these higher clusters (i.e. these will have higher number)
  # this will make their removal seem easier later hopefully
  
  
  most.likely.cluster <- c(most.likely.cluster,removed.clusters)
  most.likely.cluster <- most.likely.cluster[order(as.numeric(most.likely.cluster))]
  names.muts          <- names(most.likely.cluster)
  most.likely.cluster <- as.numeric(most.likely.cluster)
  names(most.likely.cluster) <- names.muts
  
}
if(length(most.likely.cluster)==0)
{
  
  
  solid.cluster.end <- 0
  new.cluster.start <- 1
  
  removed.clusters  <- c()    
  for (cluster in names(table(removed.muts)[order(table(removed.muts),decreasing=TRUE)]))
  {
    new.cluster <- removed.muts[removed.muts%in%cluster]
    new.cluster[new.cluster%in%cluster] <- new.cluster.start
    removed.clusters <- c(removed.clusters,new.cluster)
    new.cluster.start <- new.cluster.start +1
  }
  most.likely.cluster <- removed.clusters
  
}



#Let's add the muts to remove back into the table
v.pvals                                                        <- rep(NA,nrow(mut.table))
names(v.pvals)                                                 <- mut.table$mutation_id
v.pvals[intersect(names(v.pvals),rownames(mut.pvals))]         <- mut.pvals[intersect(names(v.pvals),rownames(mut.pvals)),1]
v.remove                                                       <- names(v.pvals)%in%muts.to.remove
names(v.remove)                                                <- names(v.pvals)
v.cluster                                                      <- rep(NA,nrow(mut.table))     
names(v.cluster)                                               <- names(v.pvals)
v.cluster[names(most.likely.cluster)]                          <- most.likely.cluster
v.size.remove                                                  <- names(v.pvals)%in%cluster.size.remove
names(v.size.remove)                                           <- names(v.pvals)
v.minCPN                                                       <- rep(NA,nrow(mut.table))
names(v.minCPN)                                                <- names(v.pvals)
v.minCPN[intersect(names(v.minCPN),rownames(mut.pvals))]       <- mut.pvals[intersect(names(v.minCPN),rownames(mut.pvals)),3]
v.regionLoss                                                      <- rep(NA,nrow(mut.table))
names(v.regionLoss)                                                <- names(v.pvals)
v.regionLoss[intersect(names(v.regionLoss),rownames(mut.pvals))]<- mut.pvals[intersect(names(v.regionLoss),rownames(mut.pvals)),4]




mut.table$cpn.remove.pval <- v.pvals
mut.table$cpn.remove      <- v.remove
mut.table$cluster         <- v.cluster
mut.table$cluster.remove  <- v.size.remove
mut.table$minCPN          <- v.minCPN
mut.table$regionLoss      <- v.regionLoss

mut.table.save.name <- paste(new.dir,patient,'.all.SNV.cpn.xls',sep="")
write.table(mut.table
            ,file=mut.table.save.name
            ,sep="\t"
            ,col.names=NA)

# let's plot these mutations ####
if(length(muts.to.remove)>1)
{
  pdf(paste(new.dir, sample, ".removedCPN.muts.pdf",sep=""),width=8,height = 8)
  clusters.to.plot <- most.likely.cluster[mut.table[mut.table$cpn.remove%in%TRUE,'mutation_id']]
  
  
  
  {
    # let's only plot a cluster if it has removed musted    
    
    print(clusters.to.plot)  
    
    lyout <- c()
    for (i in seq(1,length(regions.to.use)*2,by=2))
    {
      lyout <- rbind(lyout,rbind(c(rep(i,9),i+1)))
    }
    
    layout(lyout)
    
    for (region in regions.to.use)
    {
      
      
      region.earlyLate  <- phylo.region.list[[region]]                                             
      region.earlyLate  <- region.earlyLate[!is.na(region.earlyLate$phyloCCF),]
      region.earlyLate  <- region.earlyLate[region.earlyLate$mutation_id%in%muts.to.remove,,drop=FALSE]

      # Using seg file if exists for plotting
      region.seg.copy  <- seg.mat.copy[seg.mat.copy$SampleID%in%region,,drop=FALSE]
      
      if (!is.null(opt$input_seg_tsv)) {
        print("Using specified seg file for plotting")
        region.seg.copy     <- read.delim2(opt$input_seg_tsv)
        # If providing seg file, ensure the sample names match the sample names in input tsv
        if (!any(unique(seg.mat.copy.plot$SAMPLE) %in% unique(seg.mat.copy[,1]))){
          stop('Sample IDs do not match between input_tsv and input_seg_tsv')
        }
      } else {
        print("Using tsv data for plotting")
        region.seg.copy <- seg.mat.copy
      }

      region.seg.copy <- region.seg.copy %>% filter(SampleID %in% region)
      # ensure raw copy number columns are numeric:
      region.seg.copy$COPY_NUMBER_A <- as.numeric(region.seg.copy$COPY_NUMBER_A)
      region.seg.copy$COPY_NUMBER_B <- as.numeric(region.seg.copy$COPY_NUMBER_B)
      
      # Rename columns
      sub.mat.copy               <- region.seg.copy
      
      colnames(sub.mat.copy)[2]  <- 'Chromosome'
      colnames(sub.mat.copy)[3]  <- 'StartPosition'
      colnames(sub.mat.copy)[4]  <- 'EndPosition'
      
      #pdf(early.late.pdf)
      par(mar=c(0.5,5,0.5,0.2))
      par(lend=1)
      
      plot.simpleClusters.raw(seg.mat.patient=sub.mat.copy
                              ,most.likely.cluster = clusters.to.plot
                              ,TCGA.earlyLate=region.earlyLate
                              ,sub.clonal=1
      )
      
      mtexti(region,side = 2,off = 0.5)
      
      
      ds <- density(ifelse(as.numeric(region.earlyLate$mutCopyNum)>5,5,as.numeric(region.earlyLate$mutCopyNum)))
      ds1 <- ds
      ds1$x <- ds$y
      ds1$y <- ds$x
      par(mar=c(0.5,0,0.5,4))
      A <- hist(ifelse(as.numeric(region.earlyLate$mutCopyNum)>5,5,as.numeric(region.earlyLate$mutCopyNum)),breaks=seq(-0.25,6,by=0.1),plot=FALSE)
      plot(NULL, type = "n"
            , xlim = c(0, max(A$density))
            , ylim=c(-0.25,6)
            ,bty='n'
            ,xaxs='i'
            ,xaxt='n'
            ,yaxt='n'
            ,yaxs='i'
            ,xlab=""
            ,main=""
            ,ylab=""
      )
      rect(0, A$breaks[1:(length(A$breaks) - 1)], A$density, A$breaks[2:length(A$breaks)]
            ,border=TRUE,col="#CC6666")
      
      lines(ds1)
      
    }
  }
  
  for (cluster in unique(clusters.to.plot))
  {
    # let's only plot a cluster if it has removed musted    
    
    print(cluster)  
    
    lyout <- c()
    for (i in seq(1,length(regions.to.use)*2,by=2))
    {
      lyout <- rbind(lyout,rbind(c(rep(i,9),i+1)))
    }
    
    layout(lyout)
    
    for (region in regions.to.use)
    {
      
      
      region.earlyLate  <- phylo.region.list[[region]]                                             
      region.earlyLate  <- region.earlyLate[!is.na(region.earlyLate$phyloCCF),]
      region.earlyLate  <- region.earlyLate[region.earlyLate$mutation_id%in%muts.to.remove,,drop=FALSE]
      
      # Using seg file if exists for plotting
      region.seg.copy  <- seg.mat.copy[seg.mat.copy$SampleID%in%region,,drop=FALSE]

      if (!is.null(opt$input_seg_tsv)) {
        print("Using specified seg file for plotting")
        region.seg.copy     <- read.delim2(opt$input_seg_tsv)
        # If providing seg file, ensure the sample names match the sample names in input tsv
        if (!any(unique(seg.mat.copy.plot$SAMPLE) %in% unique(seg.mat.copy[,1]))){
          stop('Sample IDs do not match between input_tsv and input_seg_tsv')
        }
      } else {
        print("Using tsv data for plotting")
        region.seg.copy <- seg.mat.copy
      }

      region.seg.copy <- region.seg.copy %>% filter(SampleID %in% region)
      # ensure raw copy number columns are numeric:
      region.seg.copy$COPY_NUMBER_A <- as.numeric(region.seg.copy$COPY_NUMBER_A)
      region.seg.copy$COPY_NUMBER_B <- as.numeric(region.seg.copy$COPY_NUMBER_B)

      # Rename columns:
      sub.mat.copy               <- region.seg.copy
      colnames(sub.mat.copy)[2]  <- 'Chromosome'
      colnames(sub.mat.copy)[3]  <- 'StartPosition'
      colnames(sub.mat.copy)[4]  <- 'EndPosition'
      
      
      #pdf(early.late.pdf)
      par(mar=c(0.5,5,0.5,0.2))
      par(lend=1)
      
      plot.simpleClusters.raw(seg.mat.patient=sub.mat.copy
                              ,most.likely.cluster = most.likely.cluster
                              ,cluster=cluster
                              ,TCGA.earlyLate=region.earlyLate
                              ,sub.clonal=1
      )
      
      mtexti(region,side = 2,off = 0.5)
      
      
      ds <- density(ifelse(as.numeric(region.earlyLate$mutCopyNum)>5,5,as.numeric(region.earlyLate$mutCopyNum)))
      ds1 <- ds
      ds1$x <- ds$y
      ds1$y <- ds$x
      par(mar=c(0.5,0,0.5,4))
      A <- hist(ifelse(as.numeric(region.earlyLate$mutCopyNum)>5,5,as.numeric(region.earlyLate$mutCopyNum)),breaks=seq(-0.25,6,by=0.1),plot=FALSE)
      plot(NULL, type = "n"
            , xlim = c(0, max(A$density))
            , ylim=c(-0.25,6)
            ,bty='n'
            ,xaxs='i'
            ,xaxt='n'
            ,yaxt='n'
            ,yaxs='i'
            ,xlab=""
            ,main=""
            ,ylab=""
      )
      rect(0, A$breaks[1:(length(A$breaks) - 1)], A$density, A$breaks[2:length(A$breaks)]
            ,border=TRUE,col="#CC6666")
      
      lines(ds1)
      
    }
  }
  dev.off()
  
  
  # let's also write these to a table
  mut.table.remove <- mut.table[mut.table$mutation_id%in%muts.to.remove,,drop=FALSE]
  write.table(mut.table.remove
              ,file=paste(new.dir,patient,".removed.muts.txt",sep="")
              ,sep="\t"
              ,quote=FALSE
              ,col.names=NA)
  if(TRUE%in%c(mut.table.remove$driverCategory%in%driver.cat))
  {
    cat("You're removing driver muts!!")
    removed.drivers <- mut.table.remove[mut.table.remove$driverCategory%in%driver.cat,,drop=FALSE]
    write.table(removed.drivers
                ,file=paste(new.dir,patient,".removed.drivers.txt",sep="")
                ,quote=FALSE
                ,sep="\t"
                ,col.names=NA)
    cat('\n')
    cat(removed.drivers$Gene.refGene)
  }
  
}

if(ncol(mut.heatmap)>1 & !all(mut.heatmap[1] == mut.heatmap))
{
  pdf(paste(new.dir,patient,".heatmap.pdf",sep=""),height = 18,width=10)
  if(ncol(mut.heatmap)==2)
  {
    heatmap.2(mut.heatmap[order(as.numeric(rowSums(cbind(mut.heatmap[,1]+ifelse(mut.heatmap[,1]>0,2,0),mut.heatmap[,2]+ifelse(mut.heatmap[,1]>0,1,0)))),decreasing=TRUE),]
              ,Rowv = FALSE
              ,RowSideColors = ifelse(rownames(mut.heatmap)%in%muts.to.remove,'blue','yellow')
              ,col=c('gray','blue')
              ,margins=c(25,20)
              ,trace='none'
    )
              
    
    
  }
  if(ncol(mut.heatmap)>2)
  {
    heatmap.2(mut.heatmap,
              ,RowSideColors = ifelse(rownames(mut.heatmap)%in%muts.to.remove,'blue','yellow')
              ,col=c('gray','blue')
              ,margins=c(25,20)
              ,trace='none'
    )
  }
  
  dev.off()
}


no.optima = length(unique(most.likely.cluster))
max.cols = 12
require(RColorBrewer)
cols           = paste(brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
cols.opac      = paste(cols,'99',sep="")

clean.most.likely.clusters <- most.likely.cluster[most.likely.cluster%in%c(1:solid.cluster.end)]


# let's get confidence intervals for each of the mutations 
# (and also plot the confidence intervals of the driver mutations)

# # first, let's import the trace files
region.trace              <- list()
region.preClustPosterior  <- list()
region.postClustPosterior <- list()

for (region in regions.to.use) {
  clusters.to.use <- sapply(simpleClusterList, function(x) region %in% x$RegionsInCluster)
  clusters.to.use <- names(clusters.to.use)[clusters.to.use]

  for (clust in clusters.to.use) {
    pyclone.tsv <- read.table(paste0(new.dir, patient, "_cluster", clust, "/", patient, '.results.tsv')
                              , stringsAsFactors = FALSE
                              , header = TRUE)

    mutation_ids <- unlist(phylo.region.list[[region]]$mutation_id)
    if (length(grep("mutation_id", colnames(pyclone.tsv))) > 0) {
      tmp          <- intersect(mutation_ids, pyclone.tsv$mutation_id)
    } else {
      tmp          <- intersect(mutation_ids, rownames(pyclone.tsv))
    }
    phylo.region.list[[region]]$phyloCCF_PyClone.cluster[mutation_ids%in%tmp] <- most.likely.cluster[tmp]
  }
}
    

# Finally, let's put this into a megatable, and write this to an appropriate place
save(phylo.region.list,file=paste(new.dir,sample,'.PhyloRegionList.RData',sep=""))
save.image(file=paste(new.dir, sample, ".PyClone.RData",sep=""))


print("Creating human readable output")
### creating human readable output
print("Running non-original output")
tmp.phylo.region.list <- lapply(phylo.region.list, function(x) {
  tmp <- data.frame(x, stringsAsFactors = FALSE)
  rownames(tmp) <- NULL
  tmp <- tmp %>% 
    dplyr::select(mutation_id, region, Reference_Base, Alternate_Base, ref_counts, var_counts, phyloCCF, phyloCCF.0.05, phyloCCF.0.95, absolute.ccf, mutCopyNum, major_cn, minor_cn) %>%
    dplyr::mutate(mutation_id = unlist(mutation_id), ref_counts = unlist(ref_counts), var_counts = unlist(var_counts), minor_cn = unlist(minor_cn), major_cn = unlist(major_cn), region = unlist(region), Reference_Base = unlist(Reference_Base), Alternate_Base = unlist(Alternate_Base)) %>%
    dplyr::rename(SAMPLE = region, REF = Reference_Base, ALT = Alternate_Base, REF_COUNT = ref_counts, VAR_COUNT = var_counts, CCF_PHYLO = phyloCCF, CCF_OBS = absolute.ccf, MUT_COPY = mutCopyNum, COPY_NUMBER_A = major_cn, COPY_NUMBER_B = minor_cn)
  return(tmp)
})
output_tsv <- dplyr::bind_rows(tmp.phylo.region.list) 

output_tsv <- output_tsv %>% 
  dplyr::mutate(CLUSTER = most.likely.cluster[output_tsv$mutation_id],
                CLEAN = ifelse(output_tsv$mutation_id %in% names(clean.most.likely.clusters), TRUE, FALSE))

output_tsv <- output_tsv %>% 
  dplyr::mutate(CHR = sapply(strsplit(unlist(output_tsv$mutation_id), split = ":"), function(x) x[2]),
                POS = sapply(strsplit(unlist(output_tsv$mutation_id), split = ":"), function(x) x[3]),
                key = paste(paste0("chr", CHR), POS, REF, ALT, sep = ":")) %>%
  dplyr::left_join(input_tsv %>% select(key, CASE_ID, SAMPLE, DEPTH, ACF, PLOIDY), by = c("key", "SAMPLE")) %>%
  dplyr::select(CASE_ID, SAMPLE, CHR, POS, REF, ALT, REF_COUNT, VAR_COUNT, DEPTH, CLUSTER, CCF_PHYLO, CCF_OBS, MUT_COPY, COPY_NUMBER_A, COPY_NUMBER_B, ACF, PLOIDY, CLEAN, phyloCCF.0.05, phyloCCF.0.95)

write.table(output_tsv %>% select(-phyloCCF.0.05, -phyloCCF.0.95), file = paste0(new.dir, sample, ".SCoutput.FULL.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")

output_tsv_clean <- output_tsv %>% filter(CLEAN) %>% select(-CLEAN, -phyloCCF.0.05, -phyloCCF.0.95)
output_tsv_dirty <- output_tsv %>% filter(!CLEAN) %>% select(-CLEAN, -phyloCCF.0.05, -phyloCCF.0.95)

write.table(output_tsv_clean, file = paste0(new.dir, sample, ".SCoutput.CLEAN.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(output_tsv_dirty, file = paste0(new.dir, sample, ".SCoutput.DIRTY.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")


#finally, save a version of the table that is cleaned

if(opt$clean_clusters %in% TRUE)
{
  print("Cleaning clusters")
  output_tsv <- correct.clusters.from.table(output_tsv) #kg: merging clusters if sc cn correction created a new cluster
  dirty.clusters.remove <- output_tsv %>% 
    dplyr::group_by(CLUSTER) %>% 
    dplyr::filter(all(!CLEAN)) %>% 
    dplyr::pull(CLUSTER) %>% 
    unique()
  output_tsv <- output_tsv %>% 
    dplyr::mutate(CLEAN = ifelse(!CLUSTER %in% dirty.clusters.remove, TRUE, FALSE))

  ### merge clusters present in all regions at given thresholds
  clonal_cutOff        <- opt$clonal_cutOff
  propClonal_threshold <- opt$propClonal_threshold


  print("Final merging of ubiqquitous clusters")
  ### select ubiquitous mutations in clean clusters
  issue_mutations <- output_tsv %>% 
    dplyr::mutate(mutation_id = paste(CASE_ID, CHR, POS, REF, sep = ":")) %>% 
    dplyr::filter(CLEAN, mutation_id %in% ITH1muts)
  if (nrow(issue_mutations) == 0) {
    print("No additional clusters corrected")
  } else {
    clusters_to_consider <- unique(issue_mutations$CLUSTER)
    ### calculate proportion of mutations that are greater or equal to the clonal threshold in each region
    clonalProportion.df <- issue_mutations %>% 
      dplyr::group_by(CLUSTER, SAMPLE) %>%
      dplyr::mutate(nMuts = n(),
                    propClonal = sum(phyloCCF.0.95 >= clonal_cutOff) / nMuts) %>%
      dplyr::select(CLUSTER, SAMPLE, propClonal) %>%
      unique() %>% 
      dplyr::ungroup()

    ### extract the lowest proportion and filter for clusters where this is greater than the threshold
    clusters_to_change <- as.character(clonalProportion.df %>%
      dplyr::group_by(CLUSTER) %>% 
      dplyr::mutate(minPropClonal = min(propClonal)) %>%
      dplyr::select(CLUSTER, minPropClonal) %>%
      unique() %>%
      dplyr::filter(minPropClonal > propClonal_threshold) %>%
      dplyr::pull(CLUSTER))

    ### if two or more clusters are above this threshold merge all clusters into the lowest cluster ID
    if (length(clusters_to_change) == 0) {
      print("No clusters above specified thresholds")
    } else if (length(clusters_to_change) == 1) {
      print("Only single cluster above thresholds identified. Nothing to merge")
    } else {
      output_tsv <- output_tsv %>%
        dplyr::mutate(CLUSTER = ifelse(CLUSTER %in% clusters_to_change, min(clusters_to_change), CLUSTER))
    }
  }

  write.table(output_tsv %>% select(-phyloCCF.0.05, -phyloCCF.0.95), file = paste0(new.dir, sample, ".SCoutput.FULL.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
  output_tsv_clean <- output_tsv %>% filter(CLEAN) %>% select(-CLEAN, -phyloCCF.0.05, -phyloCCF.0.95)
  output_tsv_dirty <- output_tsv %>% filter(!CLEAN) %>% select(-CLEAN, -phyloCCF.0.05, -phyloCCF.0.95)
  write.table(output_tsv_clean, file = paste0(new.dir, sample, ".SCoutput.CLEAN.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(output_tsv_dirty, file = paste0(new.dir, sample, ".SCoutput.DIRTY.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")

}



if(ncol(mut.heatmap)>1)
{
  pdf(paste(new.dir,patient,".cluster.ccf.heatmap.pdf",sep=""),height = 18,width=10)
  heatmap.2(apply(heatmap.ccf,2,as.numeric),dendrogram = 'column'
            ,Rowv = FALSE
            ,RowSideColors=cols[match(most.likely.cluster[rownames(heatmap.ccf)], c(1:no.optima))]
            ,col= colorpanel(100,'white','steelblue4')
            ,margins=c(20,10)
            ,trace='none'
            ,key=FALSE
  )
  dev.off()
  
  
  pdf(paste(new.dir,patient,".cpn.ccf.heatmap.pdf",sep=""),height = 18,width=10)
  heatmap.2(apply(heatmap.ccf,2,as.numeric),dendrogram = 'column'
            ,Rowv = FALSE
            ,RowSideColors=ifelse(rownames(heatmap.ccf)%in%muts.to.remove,'red','gray')
            ,col= colorpanel(100,'white','steelblue4')
            ,margins=c(20,10)
            ,trace='none'
            ,key=FALSE
  )
  dev.off()
  
  # plot the clusters ####
  
  
  pdf(paste(new.dir,"/",patient,"_pyclone_cluster_assignment_dirty",".pdf",sep=""),height=4,width=4)
  plot.pyclone.clusters(regionList=phylo.region.list
                        ,mostLikelyClusters = most.likely.cluster
                        ,driverCat = driver.cat
                        ,ccf='absolute'
  )
  dev.off()
  
  pdf(paste(new.dir,"/",patient,"_pyclone_cluster_assignment_phylo_dirty",".pdf",sep=""),height=4,width=4)
  plot.pyclone.clusters(regionList=phylo.region.list
                        ,mostLikelyClusters = most.likely.cluster
                        ,driverCat = driver.cat
                        ,ccf='phylo'
  )
  dev.off()
  
  allMuts <- rep(1,length(most.likely.cluster))
  names(allMuts) <- names(most.likely.cluster)
  pdf(paste(new.dir,"/",patient,"_mutationCCF_all",".pdf",sep=""),height=4,width=4)
  plot.pyclone.clusters(regionList=phylo.region.list
                        ,mostLikelyClusters = allMuts
                        ,driverCat = driver.cat
                        ,ccf='phylo'
  )
  dev.off()
  
  
  
  pdf(paste(new.dir,"/",patient,"_pyclone_cluster_assignment_mutCpn_dirty",".pdf",sep=""),height=8,width=8)
  plot.pycloneMutCpn.clusters(regionList=phylo.region.list
                        ,mostLikelyClusters = most.likely.cluster
                        ,driverCat = driver.cat
                        ,ccf='phylo'
  )
  dev.off()
  
  
  pdf(paste(new.dir,"/",patient,"_pyclone_cluster_assignment_mutCpn_clean",".pdf",sep=""),height=8,width=8)
  plot.pycloneMutCpn.clusters(regionList=phylo.region.list
                              ,mostLikelyClusters = clean.most.likely.clusters
                              ,driverCat = driver.cat
                              ,ccf='phylo'
  )
  dev.off()

  pdf(paste(new.dir,"/",patient,"_pyclone_cluster_assignment_phylo_clean",".pdf",sep=""),height=4,width=4)
  plot.pyclone.clusters(regionList=phylo.region.list
                        ,mostLikelyClusters = clean.most.likely.clusters
                        ,driverCat = driver.cat
                        ,ccf='phylo'
  )
  dev.off()
  
  
  
  pdf(paste(new.dir,"/",patient,"_pyclone_cluster_assignment_phylo_clean",".pdf",sep=""),height=4,width=4)
  plot.pyclone.clusters(regionList=phylo.region.list
                        ,mostLikelyClusters = clean.most.likely.clusters
                        ,driverCat = driver.cat
                        ,ccf='phylo'
  )
  dev.off()
  
  pdf(paste(new.dir,"/",patient,"_pyclone_cluster_assignment_ccf_clean",".pdf",sep=""),height=4,width=4)
  plot.pyclone.clusters(regionList=phylo.region.list
                        ,mostLikelyClusters = clean.most.likely.clusters
                        ,driverCat = driver.cat
                        ,ccf='absolute'
  )
  dev.off()
}

# let's plot the copy number clusters
# Using seg file if exists for plotting

if (!is.null(opt$input_seg_tsv)) {
  print("Using specified seg file for plotting")
  seg.mat.copy.plot     <- read.delim2(opt$input_seg_tsv, stringsAsFactors = F)
  # If providing seg file, ensure the sample names match the sample names in input tsv
  if (!any(unique(seg.mat.copy.plot$SAMPLE) %in% unique(seg.mat.copy[,1]))){
    stop('Sample IDs do not match between input_tsv and input_seg_tsv')
  }
} else {
  print("Using tsv data for plotting")
  seg.mat.copy.plot <- seg.mat.copy
}

# ensure COPY_NUMBER_A and COPY_NUMBER_B are numeric
seg.mat.copy.plot$COPY_NUMBER_A <- as.numeric(seg.mat.copy.plot$COPY_NUMBER_A)
seg.mat.copy.plot$COPY_NUMBER_B <- as.numeric(seg.mat.copy.plot$COPY_NUMBER_B)

pdf(paste(new.dir,"/",patient,"_pyclone_cluster_assignment_copynumber_dirty",".pdf",sep=""), height = ifelse(length(regions) == 1, 5, length(regions)))
plot.region.mutCopyNum(phylo.region.list = phylo.region.list,seg.mat.copy = seg.mat.copy.plot,mostLikelyClusters = most.likely.cluster,plot.separate.clusters = TRUE)
dev.off()

pdf(paste(new.dir,"/",patient,"_pyclone_cluster_assignment_copynumber_clean",".pdf",sep=""), height = ifelse(length(regions) == 1, 5, length(regions)))
plot.region.mutCopyNum(phylo.region.list = phylo.region.list,seg.mat.copy = seg.mat.copy.plot,mostLikelyClusters = clean.most.likely.clusters,plot.separate.clusters = TRUE)
dev.off()


save.image(file=paste(new.dir, sample, ".PyClone.RData",sep=""))


print(sessionInfo());


pyclone.out.complete <- paste(new.dir,patient,".simpleClustering.complete",sep='');
system(paste("touch ",pyclone.out.complete,sep=""));


