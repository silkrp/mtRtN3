/*
 * Run RtN on chrM bam file to remove NUMT reads
 */
process RTN {
    
    label 'bigmem'
    container 'rpsilk/ge-mtdna'

    input:
    tuple val(meta), path(bam), path(index)
    path(hmtdb)
    path(numts)

    output:
    tuple val(meta), path("${meta.sample}_chrm.rtn.bam"), emit: numt_filt_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    def group = task.ext.prefix ?: "${meta.group}"

    """
    rtn -p -h humans.fa -n Calabrese_Dayama_Smart_Numts.fa -c ${params.contig} -b $bam
    """

    stub:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    def group = task.ext.prefix ?: "${meta.group}"
    
    """
    touch ${sample}_chrm.rtn.bam
    """
}
