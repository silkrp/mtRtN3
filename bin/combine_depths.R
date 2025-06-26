#!/usr/bin/env Rscript

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine sequencing depths for individual mtDNA base-pairs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(data.table)

args <- commandArgs(trailingOnly = T)

out <- data.frame(chr = "chrM", pos = seq(1,16569))
for(i in args[-length(args)]){
  
  sample <- strsplit(basename(i), "[.]")[[1]][1]
  
  if(args[length(args)] == "true"){
    
    group = strsplit(basename(i), "[.]")[[1]][2]
    sample = paste0(sample, ".", group)
    
  }
  
  single <- fread(i, header = F)
  colnames(single) <- c("chr", "pos", sample)
  single$chr <- "chrM"
  
  out <- merge(out, single, by = c("chr", "pos"))
  
}

fwrite(out, "depths.txt.gz", sep = "\t", quote = F, row.names = F)