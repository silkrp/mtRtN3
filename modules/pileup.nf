/*
 * Generate pileup file using samtools
 */
process PILEUP {
    
    label 'process_small'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)
    path(fasta)

    output:
    tuple val(meta), path("${meta.sample}.${meta.group}.pileup"), emit: mpileup_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    def group = task.ext.prefix ?: "${meta.group}"

    """
    samtools sort $bam | \
        samtools mpileup -d 0 -f $fasta - > ${sample}.${group}.pileup
    """

    stub:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    def group = task.ext.prefix ?: "${meta.group}"
    
    """
    touch ${sample}.${group}.pileup
    """
}