#!/usr/bin/env Rscript

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate mtDNA copy number for each sample
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(data.table)
library(dplyr)

# read in command line arguments 
args <- commandArgs(trailingOnly = TRUE)
#args <- c("/Users/rsilk/Desktop/TCGA-ZN-A9VW.tumour.mosdepth.summary.txt", "TCGA-ZN-A9VW", "false", "false", "tumour", "chrM")

depths <- fread(args[1])
sample <- args[2]
purity <- args[3]
ploidy <- args[4]
type <- args[5]
contig <- args[6]

auto <- filter(depths, chrom == "total")
chrm <- filter(depths, chrom == contig)

if (purity == "false" & ploidy == "false") {
  
  cov_chrm <- chrm$mean
  cov_auto <- auto$mean
  
  copy_number_simple <- round((cov_chrm / cov_auto), 2)
  copy_number <- round((cov_chrm / cov_auto) * (0 * 2 + (1 - 0) * 2),2)
  
  out <- data.frame(sampleID = sample, 
                    type = type,
                    mtDNA_avg_depth = cov_chrm, 
                    auto_avg_depth = cov_auto, 
                    copy_number_simple = copy_number_simple, 
                    copy_number = copy_number)
  
  out_name <- paste0(sample, ".", type, ".cn.txt")
  fwrite(out, out_name, sep = "\t", quote = F, row.names = F)
  
} else if (purity == "false" & ploidy != "false") {
  
  ploidy <- fread(ploidy)
  
  cov_chrm <- chrm$mean
  cov_auto <- auto$mean
  
  if(type == "tumour") {
    
    ploidy_est <- dplyr::filter(ploidy, sampleID == sample)
    print(class(ploidy_est))
    ploidy_est <- ploidy_est$ploidy
    
  } else {
    
    ploidy_est <- 2
    
  }

  copy_number_simple <- round((cov_chrm / cov_auto), 2)
  copy_number <- round((cov_chrm / cov_auto) * (1 * ploidy_est + (1 - 1) * 2),2)
  
  out <- data.frame(sampleID = sample, 
                    type = type,
                    mtDNA_avg_depth = cov_chrm, 
                    auto_avg_depth = cov_auto, 
                    copy_number_simple = copy_number_simple, 
                    copy_number = copy_number)
  
  out_name <- paste0(sample, ".", type, ".cn.txt")
  fwrite(out, out_name, sep = "\t", quote = F, row.names = F)
  
} else if (purity != "false" & ploidy == "false") {
  
  purity <- fread(purity)
  
  cov_chrm <- chrm$mean
  cov_auto <- auto$mean
  
  purity_est <- filter(purity, sampleID == sample & group == type)
  
  if(type == "tumour"){
    purity_est <- purity_est$purity
  } else {
    purity_est <- 1 - purity_est$purity
  }
  
  copy_number_simple <- round((cov_chrm / cov_auto), 2)
  copy_number <- round((cov_chrm / cov_auto) * (purity_est * 2 + (1 - purity_est) * 2),2)
  
  out <- data.frame(sampleID = sample, 
                    type = type,
                    mtDNA_avg_depth = cov_chrm, 
                    auto_avg_depth = cov_auto, 
                    copy_number_simple = copy_number_simple, 
                    copy_number = copy_number)
  
  out_name <- paste0(sample, ".", type, ".cn.txt")
  fwrite(out, out_name, sep = "\t", quote = F, row.names = F)
  
} else {
  
  ploidy <- fread(ploidy)
  purity <- fread(purity)

  cov_chrm <- chrm$mean
  cov_auto <- auto$mean
  
  if(type == "tumour") {
    
    ploidy_est <- filter(ploidy, sampleID == sample)
    print(ploidy_est)
    ploidy_est <- ploidy_est$ploidy
    
  } else {
    
    ploidy_est <- 2
    
  }
  
  purity_est <- filter(purity, sampleID == sample & group == type)
  
  if(type == "tumour"){
    purity_est <- purity_est$purity
  } else {
    purity_est <- 1 - purity_est$purity
  }

  copy_number_simple <- round((cov_chrm / cov_auto), 2)
  copy_number <- round((cov_chrm / cov_auto) * (purity_est * ploidy_est + (1 - purity_est) * 2),2)
  
  out <- data.frame(sampleID = sample,
                    type = type,
                    mtDNA_avg_depth = cov_chrm, 
                    auto_avg_depth = cov_auto, 
                    copy_number_simple = copy_number_simple, 
                    copy_number = copy_number)
  
  out_name <- paste0(sample, ".", type, ".cn.txt")
  fwrite(out, out_name, sep = "\t", quote = F, row.names = F)
  
}

      