/*
 * Get read depths for each base-pair in the mtDNA 
 */
process MTDNA_DEPTHS {
    
    label 'process_small'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(index)

    output:
    path("${meta.sample}.${meta.group}.depths.txt"), emit: depths_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    def group = task.ext.prefix ?: "${meta.group}"

    """
    samtools depth -d 1000000 -r ${params.contig} -aa $bam > ${sample}.${group}.depths.txt
    """

    stub:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    def group = task.ext.prefix ?: "${meta.group}"
    
    """
    touch ${sample}.${group}.depths.txt
    """
}
