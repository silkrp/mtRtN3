/*
 * Get the average sequencing depth across the mtDNA and the nuclear genome
 */
process COVERAGE {
    
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.8--hd299d5a_0' :
        'biocontainers/mosdepth:0.3.8--hd299d5a_0'}"

    input:
    tuple val(meta), path(bam), path(index)
    path(cram_fasta)

    output:
    tuple val(meta), path("${meta.sample}.${meta.group}.mosdepth.summary.txt"), emit: coverage_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    def group = task.ext.prefix ?: "${meta.group}"
    //def fasta = cram_fasta ? "--fasta ${cram_fasta}" : ""
    def fasta = cram_fasta.name != 'NO_FILE' ? "--fasta ${cram_fasta}" : ''

    """
    mosdepth $fasta ${sample}.${group} $bam
    """

    stub:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    def group = task.ext.prefix ?: "${meta.group}"
    
    """
    touch ${sample}.${group}.mosdepth.summary.txt
    """
}
