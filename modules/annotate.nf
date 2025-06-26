/*
 * Combine VEP annotations and mtDNA variant calls
 */
process ANNOTATE {
    
    label 'small_process'
    container 'rpsilk/rdtds'
    publishDir "${params.outdir}/results/" , mode:'copy'

    input:
    tuple path(snvs), path(indels)
    tuple path(snvs_vep), path(indels_vep)
    path(gnomad_vcf)

    output:
    path("allcalls.txt.gz"), emit: final_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    annotate.R $snvs $snvs_vep $indels $indels_vep $gnomad_vcf $params.matched
    """

    stub:    
    """
    touch allcalls.txt.gz
    """
}

