# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Filter mtRtN3 Variant Calls
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(data.table)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(MutationalPatterns)
library(ggplot2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define functions to filter variants and check mtDNA mutational signature
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

filter_varscan <- function(raw_varscan_calls, min_supporting_reads = 10, min_vaf = 1, remove_biallelic = FALSE, matched = TRUE){
  
  if(matched){
    
    filtered <- raw_varscan_calls %>%
      mutate(tumor_var_freq = as.numeric(sub("%", "", tumor_var_freq)),
             normal_var_freq = as.numeric(sub("%", "", normal_var_freq))) %>%
      
      # remove homopolymer variants
      filter((position <= 66 | position >= 71) & 
               (position <= 300 | position >= 316) & 
               (position <= 513 | position >= 525) &
               (position <= 3106 | position >= 3107) & 
               (position <= 12418 | position >= 12425) &
               (position <= 16182 | position >= 16194)) %>%
      
      # remove variants with low VAF and low supporting reads
      filter((somatic_status == "Somatic" & tumor_reads2 >= min_supporting_reads & tumor_var_freq >= min_vaf) |
               (somatic_status == "Germline" & normal_reads2 >= min_supporting_reads & normal_var_freq >= min_vaf))
    
    # calculate strand bias
    filtered$strand_bias <- apply(filtered, 1, function(x){ ifelse(x[13] == "Somatic",
                                                                   fisher.test(matrix(as.numeric(x[16:19]), ncol = 2, byrow = T))$p.value, 
                                                                   fisher.test(matrix(as.numeric(x[20:23]), ncol = 2, byrow = T))$p.value) })
    filtered$phred_score <- -10 * log10(filtered$strand_bias)
    filtered <- filter(filtered, phred_score < 60)
    
    # remove biallelic positions
    if(remove_biallelic){
      
      biallelic_positions <- filtered %>% group_by(position, sample) %>% tally() %>% filter(n > 1)
      filtered <- filtered %>% filter(!(position %in% biallelic_positions$position) & !(sample %in% biallelic_positions$sample))
      
    }
    
    return(filtered)

  } else {
    
    filtered <- raw_varscan_calls %>%
      mutate(Freq = as.numeric(sub("%", "", Freq))) %>%
      
      # remove homopolymer variants
      filter((Position <= 66 | Position >= 71) & 
               (Position <= 300 | Position >= 316) & 
               (Position <= 513 | Position >= 525) &
               (Position <= 3106 | Position >= 3107) & 
               (Position <= 12418 | Position >= 12425) &
               (Position <= 16182 | Position >= 16194)) %>%
      
      # remove variants with low VAF and low supporting reads
      filter(Reads2 >= min_supporting_reads & Freq >= min_vaf)
    
    # calculate strand bias
    filtered$strand_bias <- apply(filtered, 1, function(x){ fisher.test(matrix(as.numeric(x[12:15]), ncol = 2, byrow = T))$p.value })
    filtered$phred_score <- -10 * log10(filtered$strand_bias)
    filtered <- filter(filtered, phred_score < 60)
    
    # remove biallelic positions
    if(remove_biallelic){
      
      biallelic_positions <- filtered %>% group_by(Position, sample) %>% tally() %>% filter(n > 1)
      filtered <- filtered %>% filter(!(Position %in% biallelic_positions$position) & !(sample %in% biallelic_positions$sample))
      
    }
    
    return(filtered)
    
  }
}

plot_mut_sig = function(filtered_variants, matched = TRUE){
  
  if(matched){
    
    muts <- filtered_variants %>%
      dplyr::filter(mut_type == "SNV" & somatic_status == "Somatic") %>%
      dplyr::mutate(strand = ifelse(ref == "C", "-",
                                    ifelse(ref == "G", "+",
                                           ifelse(ref == "T", "-",
                                                  ifelse(ref == "A", "+", "OTHER")))))
    
    light <- filter(muts, strand == "-")
    heavy <- filter(muts, strand == "+")
    
    light_gr <- with(light, GRanges(seqnames = Rle(chrom), ranges = IRanges(start = position, end = position), REF = ref, ALT = var, strand = strand))
    heavy_gr <- with(heavy, GRanges(seqnames = Rle(chrom), ranges = IRanges(start = position, end = position), REF = ref, ALT = var, strand = strand))
    
    GenomeInfoDb::genome(light_gr) = 'hg38'
    GenomeInfoDb::genome(heavy_gr) = 'hg38'
    
    gL <- GRangesList("light" = light_gr,
                      "heavy" = heavy_gr)
    
    mut_matrix <- mut_matrix(vcf_list = gL, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
    
    mut_matrix <- as.data.frame(mut_matrix)
    mut_matrix$id <- rownames(mut_matrix)
    
    mut_matrix <- reshape2::melt(mut_matrix, id = "id")
    
    mut_matrix$groups <- substring(mut_matrix$id, 3, 5)
    mut_matrix$xlab <- paste0(substring(mut_matrix$id, 1, 1), "-",substring(mut_matrix$id, 7, 7))
    
    plot <- ggplot(mut_matrix) +
      geom_bar(aes(x = xlab, y = value, fill = factor(variable)), stat = "identity") + 
      facet_wrap(~groups, nrow = 1, ncol = 6) +
      scale_fill_brewer(palette = "Set1", name = "Strand", labels = c("Light", "Heavy")) +
      theme_gray() +
      labs(x = "", y = "Mutations\n") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
            panel.spacing = unit(0.1, "lines"),
            panel.grid.major.x = element_blank(),
            legend.key.size = unit(0.4, 'cm'),
            legend.title = element_text(size=10),
            legend.text = element_text(size=9)) 
    
    return(plot)
    
  } else {
    
    muts <- filtered_variants %>%
      dplyr::filter(mut_type == "SNV") %>%
      dplyr::mutate(strand = ifelse(Ref == "C", "-",
                                    ifelse(Ref == "G", "+",
                                           ifelse(Ref == "T", "-",
                                                  ifelse(Ref == "A", "+", "OTHER")))))
    
    light <- filter(muts, strand == "-")
    heavy <- filter(muts, strand == "+")
    
    light_gr <- with(light, GRanges(seqnames = Rle(Chrom), ranges = IRanges(start = Position, end = Position), REF = Ref, ALT = Var, strand = strand))
    heavy_gr <- with(heavy, GRanges(seqnames = Rle(Chrom), ranges = IRanges(start = Position, end = Position), REF = Ref, ALT = Var, strand = strand))
    
    GenomeInfoDb::genome(light_gr) = 'hg38'
    GenomeInfoDb::genome(heavy_gr) = 'hg38'
    
    gL <- GRangesList("light" = light_gr,
                      "heavy" = heavy_gr)
    
    mut_matrix <- mut_matrix(vcf_list = gL, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
    
    mut_matrix <- as.data.frame(mut_matrix)
    mut_matrix$id <- rownames(mut_matrix)
    
    mut_matrix <- reshape2::melt(mut_matrix, id = "id")
    
    mut_matrix$groups <- substring(mut_matrix$id, 3, 5)
    mut_matrix$xlab <- paste0(substring(mut_matrix$id, 1, 1), "-",substring(mut_matrix$id, 7, 7))
    
    plot <- ggplot(mut_matrix) +
      geom_bar(aes(x = xlab, y = value, fill = factor(variable)), stat = "identity") + 
      facet_wrap(~groups, nrow = 1, ncol = 6) +
      scale_fill_brewer(palette = "Set1", name = "Strand", labels = c("Light", "Heavy")) +
      theme_gray() +
      labs(x = "", y = "Mutations\n") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 5),
            panel.spacing = unit(0.1, "lines"),
            panel.grid.major.x = element_blank(),
            legend.key.size = unit(0.4, 'cm'),
            legend.title = element_text(size=10),
            legend.text = element_text(size=9))
    
    return(plot)
    
  }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Apply filtering and plot mutational signature to raw variant calls
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

filtered <- fread("/Users/rsilk/Documents/twist/results/allcalls.txt.gz") %>%
  filter_varscan(min_supporting_reads = 10, min_vaf = 1, remove_biallelic = FALSE, matched = FALSE)

nrow(filtered)
table(filtered$mut_type)

plot_mut_sig(filtered, matched = FALSE)

fwrite(filtered, "/Users/rsilk/Documents/twist/results/allcallsfilt.txt.gz", sep = "\t", quote = F, row.names = F)
