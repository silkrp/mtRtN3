/*
 * Run VarScan2 on pileup file of mtDNA reads
 */
process VARSCAN2 {
    
    label 'process_small'
    container 'rpsilk/ge-mtdna'

    input:
    tuple val(meta), path(t_pileup), path(n_pileup)

    output:
    path("${meta.sample}.snp"), emit: snvs_ch
    path("${meta.sample}.indel"), emit: indels_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"

    if ( "${params.matched}" == "true" ){
        """
        varscan somatic ${n_pileup} ${t_pileup} ${sample} \
            --strand-filter 1 --min-avg-qual 30 --min-coverage 2 --min-reads 2 --min-var-freq 0 \
            --output-snp output.snp --output-indel output.indel
        
        awk '{print \$0 "\t" (NR==1 ? "sample" : "${sample}")}' output.snp > ${sample}.snp
        awk '{print \$0 "\t" (NR==1 ? "sample" : "${sample}")}' output.indel > ${sample}.indel
        """
    } else if ( "${params.matched}" == "false" ) {
        """
        varscan mpileup2snp ${t_pileup} \
            --strand-filter 1 --min-avg-qual 30 \
            --min-coverage 2 --min-reads2 2 --min-var-freq 0 > output.snp

        varscan mpileup2indel ${t_pileup} \
            --strand-filter 1 --min-avg-qual 30 \
            --min-coverage 2 --min-reads2 2 --min-var-freq 0 > output.indel

        awk '{print \$0 "\t" (NR==1 ? "sample" : "${sample}")}' output.snp > ${sample}.snp
        awk '{print \$0 "\t" (NR==1 ? "sample" : "${sample}")}' output.indel > ${sample}.indel
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def sample = task.ext.prefix ?: "${meta.sample}"
    
    """
    touch ${sample}.snp 
    touch ${sample}.indel
    """
}


//--normal-purity ${n_purity} --tumor-purity ${t_purity}