//
// Prepare reference genome files for mtDNA and NUMTs 
//
workflow PREPARE_GENOME {
    
    take:
    fasta                //      file: /path/to/genome.fasta
    hmtdb                //      file: /path/to/hmtdb 
    numts                //      file: /path/to/numts
    cram_fasta           //      file: /path/to/cram.reference
    vep_cache
    vep_species
    vep_genome
    vep_cache_version
    gnomad_vcf

    main:
    
    ch_fasta = Channel.value(file(fasta))
    hmtdb_ch = Channel.value(file(hmtdb))
    numts_ch = Channel.value(file(numts))

    //cram_fasta_ch = cram_fasta ? Channel.value(file(cram_fasta)) : Channel.empty()
    cram_fasta_ch = file(cram_fasta, checkIfExists:true)
    
    gnomad_ch = Channel.value(file(gnomad_vcf))

    def vep_annotation_cache_key = "${vep_cache_version}_${vep_genome}/"
    ensemblvep_cache = Channel.fromPath(file("${vep_cache}"), checkIfExists: true).collect()

    emit:
    fasta = ch_fasta                  // channel: path(genome.fasta)
    hmtdb = hmtdb_ch                  // channel: path(rtn.hmtdb)
    numts = numts_ch                  // channel: path(rtn.numts)
    cram_fasta = cram_fasta_ch        // channel: path(cram.fasta)
    vep_cache = ensemblvep_cache
    gnomad_vcf = gnomad_ch

}
