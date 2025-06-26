/*
 * Get read depths for each base-pair in the mtDNA 
 */
process COMBINE_DEPTHS {
    
    label 'process_small'
    container 'rpsilk/rdtds'
    publishDir "${params.outdir}/results/" , mode:'copy'

    input:
    path(depths)

    output:
    path("depths.txt.gz")

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    combine_depths.R $depths ${params.matched}
    """

    stub:
    """
    touch depths.txt.gz
    """
}
