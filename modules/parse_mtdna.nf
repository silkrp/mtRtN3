/*
 * Extract only the mitochondria reads from the bam files
 */
process PARSE_MTDNA {
    
    label 'process_small'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(index)
    path(cram_fasta)

    output:
    tuple val(meta), path("${meta.sample}_chrm.bam"), path("${meta.sample}_chrm.bam.bai"), emit: chrm_bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    //def reference = cram_fasta ? "--reference ${cram_fasta}" : ""
    def reference = cram_fasta.name != 'NO_FILE' ? "--reference ${cram_fasta}" : ''

    """
    samtools view -h -b $reference $bam $params.contig > ${sample}_chrm.bam
    samtools index ${sample}_chrm.bam
    """

    stub:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    
    """
    touch ${sample}_chrm.bam
    touch ${sample}_chrm.bam.bai
    """
}
