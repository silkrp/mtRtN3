/*
 * Calculate the copy number for each sample
 */
process CALCULATE_CN {
    
    label 'process_small'
    container 'rpsilk/rdtds'

    input:
    tuple val(meta), path(depths)

    output:
    path("${meta.sample}.${meta.group}.cn.txt"), emit: cn_ch
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    def group = task.ext.prefix ?: "${meta.group}"

    """
    calculate_cn.R $depths $sample "false" "false" $group $params.contig
    """

    stub:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    def group = task.ext.prefix ?: "${meta.group}"
    
    """
    touch ${sample}.${group}.cn.txt
    """
}
