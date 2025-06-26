#!/usr/bin/env Rscript

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Combine VEP annotations and mtDNA variant calls
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(data.table)
library(dplyr)
library(vcfR)
library(tidyr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
#args <- c("snvs", "snvs_vep", "indels", "indels_vep", "gnomad_vcf", "matched")

## If the samples are matched run this
## else run unmatched version
if(args[6] == "true"){
  
  snvs <- fread(args[1])
  snvs_vep <- fread(args[2], skip = 107) %>%
    aggregate(data = ., .~`#Uploaded_variation`, FUN = paste, collapse = ";")
  
  indels <- fread(args[3])
  indels_vep <- fread(args[4], skip = 107) %>%
    aggregate(data = ., .~`#Uploaded_variation`, FUN = paste, collapse = ";")
  
  ### Merge VEP annotations and filtered calls
  snvs_annotated <- snvs %>% 
    mutate(`#Uploaded_variation` = paste0(chrom, "_", position, "_", ref, "/", var)) %>%
    left_join(., snvs_vep, by = "#Uploaded_variation") %>%
    dplyr::select(-`#Uploaded_variation`) %>%
    distinct()
  
  indels_annotated <- indels %>%
    mutate(`#Uploaded_variation` = paste0(position, "_", var)) %>%
    left_join(., indels_vep, by = "#Uploaded_variation") %>%
    dplyr::select(-`#Uploaded_variation`) %>%
    distinct()
  
  annotated <- rbind(snvs_annotated, indels_annotated)
  
  ### Add in additional identifiers 
  annotated <- annotated %>% 
    mutate(tag = paste0(position, "_", ref, var),
           full_tag = paste0(sample, "_", position, "_", ref, var),
           mut_type = ifelse(grepl("-", var, fixed = T), "DEL", ifelse(grepl("+", var, fixed = T), "INS", "SNV")))
  
  ### Annotate haplotype defining variants 
  gnomad <- read.vcfR(args[5])
  
  hap_defining <- gnomad@fix %>%
    as.data.frame() %>%
    filter(grepl("hap_defining_variant", INFO)) %>%
    mutate(POS = as.numeric(POS),
           tag = paste0(POS, "_", REF, ALT))
  
  annotated <- annotated %>%
    mutate(vcf_ref = ifelse(grepl("+", var, fixed = TRUE), ref, 
                            ifelse(grepl("-", var, fixed = TRUE), paste0(ref, sub(x = var, pattern = "-", "", fixed = TRUE)), ref)),
           vcf_alt = ifelse(grepl("-", var, fixed = TRUE), ref, 
                            ifelse(grepl("+", var, fixed = TRUE), paste0(ref, sub(x = var, pattern = "+", "", fixed = TRUE)), var)),
           vcf_tag = paste0(position, "_", vcf_ref, vcf_alt),
           haplotype_defining = ifelse(vcf_tag %in% hap_defining$tag, TRUE, FALSE))
  
  ### Write out the final annotated variant calls
  fwrite(annotated, "allcalls.txt.gz", sep = "\t", quote = F, row.names = F)
  
} else {
  
  snvs <- fread(args[1]) %>%
    subset(., select = which(!duplicated(names(.)))) %>%
    separate(`Cons:Cov:Reads1:Reads2:Freq:P-value`, into = c("Cons", "Cov", "Reads1", "Reads2", "Freq", "P_value"), sep = ":") %>%
    separate(`StrandFilter:R1+:R1-:R2+:R2-:pval`, into = c("StrandFilter","R1+","R1-","R2+","R2-","pval","pval2"), sep = ":") 
  
  snvs_vep <- fread(args[2], skip = 107) %>%
    aggregate(data = ., .~`#Uploaded_variation`, FUN = paste, collapse = ";")
  
  indels <- fread(args[3]) %>%
    subset(., select = which(!duplicated(names(.)))) %>%
    separate(`Cons:Cov:Reads1:Reads2:Freq:P-value`, into = c("Cons", "Cov", "Reads1", "Reads2", "Freq", "P_value"), sep = ":") %>%
    separate(`StrandFilter:R1+:R1-:R2+:R2-:pval`, into = c("StrandFilter","R1+","R1-","R2+","R2-","pval","pval2"), sep = ":") 
  
  indels_vep <- fread(args[4], skip = 107) %>%
    aggregate(data = ., .~`#Uploaded_variation`, FUN = paste, collapse = ";")
  
  ### Merge VEP annotations and filtered calls
  snvs_annotated <- snvs %>% 
    mutate(`#Uploaded_variation` = paste0(Chrom, "_", Position, "_", Ref, "/", Var)) %>%
    left_join(., snvs_vep, by = "#Uploaded_variation") %>%
    dplyr::select(-`#Uploaded_variation`) %>%
    distinct()
  
  indels_annotated <- indels %>%
    mutate(`#Uploaded_variation` = paste0(Position, "_", Var)) %>%
    left_join(., indels_vep, by = "#Uploaded_variation") %>%
    dplyr::select(-`#Uploaded_variation`) %>%
    distinct()
  
  annotated <- rbind(snvs_annotated, indels_annotated)
  
  ### Add in additional identifiers 
  annotated <- annotated %>% 
    mutate(tag = paste0(Position, "_", Ref, Var),
           full_tag = paste0(sample, "_", Position, "_", Ref, Var),
           mut_type = ifelse(grepl("-", Var, fixed = T), "DEL", ifelse(grepl("+", Var, fixed = T), "INS", "SNV")))
  
  ### Annotate haplotype defining variants 
  gnomad <- read.vcfR(args[5])
  
  hap_defining <- gnomad@fix %>%
    as.data.frame() %>%
    filter(grepl("hap_defining_variant", INFO)) %>%
    mutate(POS = as.numeric(POS),
           tag = paste0(POS, "_", REF, ALT))
  
  annotated <- annotated %>%
    mutate(vcf_ref = ifelse(grepl("+", Var, fixed = TRUE), Ref, 
                            ifelse(grepl("-", Var, fixed = TRUE), paste0(Ref, sub(x = Var, pattern = "-", "", fixed = TRUE)), Ref)),
           vcf_alt = ifelse(grepl("-", Var, fixed = TRUE), Ref, 
                            ifelse(grepl("+", Var, fixed = TRUE), paste0(Ref, sub(x = Var, pattern = "+", "", fixed = TRUE)), Var)),
           vcf_tag = paste0(Position, "_", vcf_ref, vcf_alt),
           haplotype_defining = ifelse(vcf_tag %in% hap_defining$tag, TRUE, FALSE))
  
  ### Write out the final annotated variant calls
  fwrite(annotated, "allcalls.txt.gz", sep = "\t", quote = F, row.names = F)
  
}


