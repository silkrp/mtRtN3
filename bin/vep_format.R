#!/usr/bin/env Rscript

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Format VarScan2 calls into vcf format for VEP
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(data.table)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("/Users/rsilk/Desktop/combinedSNPs.txt", "/Users/rsilk/Desktop/combinedIndels.txt", "true")

if(args[3] == "true"){
  
  SNVs <- fread(args[1])
  INDELs <- fread(args[2])
  
  vep_snvs_input <- SNVs %>%
    mutate(change = paste0(ref, "/", var), stop = position, num = 1) %>%
    dplyr::select(chrom, position, stop, change) %>%
    arrange(position) %>%
    distinct()
  
  vep_indels_input <- INDELs %>%
    mutate(vcf_ref = ifelse(grepl("+", var, fixed = TRUE), ref, 
                            ifelse(grepl("-", var, fixed = TRUE), paste0(ref, sub(x = var, pattern = "-", "", fixed = TRUE)), ref)),
           vcf_alt = ifelse(grepl("-", var, fixed = TRUE), ref, 
                            ifelse(grepl("+", var, fixed = TRUE), paste0(ref, sub(x = var, pattern = "+", "", fixed = TRUE)), var)),
           identifier = paste0(position, "_", var),
           point2 = ".", point3 = ".", point4 = ".") %>%
    dplyr::select(chrom, position, identifier, vcf_ref, vcf_alt, point2, point3, point4) %>% 
    arrange(position) %>%
    distinct()
  
  fwrite(vep_snvs_input, "snvs.input", sep = "\t", quote = F, row.names = F, col.names = F)
  fwrite(vep_indels_input, "indels.input", sep = "\t", quote = F, row.names = F, col.names = F)

} else {
  
  SNVs <- fread(args[1]) %>%
    subset(., select = which(!duplicated(names(.)))) %>%
    separate(`Cons:Cov:Reads1:Reads2:Freq:P-value`, into = c("Cons", "Cov", "Reads1", "Reads2", "Freq", "P_value"), sep = ":") %>%
    separate(`StrandFilter:R1+:R1-:R2+:R2-:pval`, into = c("StrandFilter","R1+","R1-","R2+","R2-","pval","pval2"), sep = ":") 
    
  INDELs <- fread(args[2]) %>%
    subset(., select = which(!duplicated(names(.)))) %>%
    separate(`Cons:Cov:Reads1:Reads2:Freq:P-value`, into = c("Cons", "Cov", "Reads1", "Reads2", "Freq", "P_value"), sep = ":") %>%
    separate(`StrandFilter:R1+:R1-:R2+:R2-:pval`, into = c("StrandFilter","R1+","R1-","R2+","R2-","pval","pval2"), sep = ":") 
  
  vep_snvs_input <- SNVs %>%
    mutate(change = paste0(Ref, "/", Var), stop = Position, num = 1) %>%
    dplyr::select(Chrom, Position, stop, change) %>%
    arrange(Position) %>%
    distinct()
  
  vep_indels_input <- INDELs %>%
    mutate(vcf_ref = ifelse(grepl("+", Var, fixed = TRUE), Ref, 
                            ifelse(grepl("-", Var, fixed = TRUE), paste0(Ref, sub(x = Var, pattern = "-", "", fixed = TRUE)), Ref)),
           vcf_alt = ifelse(grepl("-", Var, fixed = TRUE), Ref, 
                            ifelse(grepl("+", Var, fixed = TRUE), paste0(Ref, sub(x = Var, pattern = "+", "", fixed = TRUE)), Var)),
           identifier = paste0(Position, "_", Var),
           point2 = ".", point3 = ".", point4 = ".") %>%
    dplyr::select(Chrom, Position, identifier, vcf_ref, vcf_alt, point2, point3, point4) %>% 
    arrange(Position) %>%
    distinct()
  
  fwrite(vep_snvs_input, "snvs.input", sep = "\t", quote = F, row.names = F, col.names = F)
  fwrite(vep_indels_input, "indels.input", sep = "\t", quote = F, row.names = F, col.names = F)
  
}
