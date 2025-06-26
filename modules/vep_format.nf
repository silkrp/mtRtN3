/*
 * Format variant calls for VEP
 */
process VEP_FORMAT {
    
    label 'small_process'
    container 'rpsilk/rdtds'

    input:
    tuple path(snvs), path(indels)

    output:
    tuple path("snvs.input"), path("indels.input"), emit: vep_formatted_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    vep_format.R $snvs $indels ${params.matched}
    """

    stub:    
    """
    touch snvs.input
    touch indels.input
    """
}