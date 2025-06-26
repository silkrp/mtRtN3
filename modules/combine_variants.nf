/*
 * Merge all variant calls
 */
process COMBINE_VARIANTS {
    
    label 'process_small'
    container 'rpsilk/rdtds'

    input:
    path(indels)
    path(snps)

    output:
    tuple path("combinedSNPs.txt"), path("combinedIndels.txt"), emit: combined_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'FNR==1 && NR!=1 {next} {print}' $snps > combinedSNPs.txt
    awk 'FNR==1 && NR!=1 {next} {print}' $indels > combinedIndels.txt
    """

    stub:
    """
    touch combinedSNPs.txt
    touch combinedIndels.txt
    """
}