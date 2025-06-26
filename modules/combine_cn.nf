/*
 * Combine all mtDNA copy-number estimates
 */
process COMBINE_CN {
    
    label 'process_small'
    container 'rpsilk/rdtds'
    publishDir "${params.outdir}/results/" , mode:'copy'

    input:
    path(cn_estimates)

    output:
    path("copies.txt.gz")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk 'FNR==1 && NR!=1 {next} {print}' $cn_estimates > copies.txt
    gzip copies.txt
    """

    stub:    
    """
    touch copies.txt.gz
    """
}
