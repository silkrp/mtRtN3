/*
 * Get variant annotation using VEP
 */
process VEP {
    
    label 'bigmem'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:113.0--pl5321h2a3209d_0' :
        'biocontainers/ensembl-vep:113.0--pl5321h2a3209d_0' }"

    input:
    tuple path(snvs), path(indels)
    path fasta
    path vep_cache

    output:
    tuple path("snvs.vep.output"), path("indels.vep.output"), emit: annotated_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def dir_cache = "\${PWD}/${vep_cache}"

    """
    vep -i $snvs -o snvs.vep.output --tab --fasta $fasta --assembly GRCh38 --species homo_sapiens \
        --cache --cache_version 113 --dir_cache $dir_cache --everything --distance 0 

    vep -i $indels -o indels.vep.output --tab --fasta $fasta --assembly GRCh38 --species homo_sapiens \
        --cache --cache_version 113 --dir_cache $dir_cache --everything --distance 0 
    """

    stub:    
    """
    touch snvs.vep.output
    touch indels.vep.output
    """
}







