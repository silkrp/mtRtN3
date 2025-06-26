//
// Run variant calling pipeline
//

include { RTN } from '../modules/rtn'
include { PILEUP } from '../modules/pileup'
include { VARSCAN2 } from '../modules/varscan2'
include { COMBINE_VARIANTS } from '../modules/combine_variants'
include { VEP_FORMAT } from '../modules/vep_format'
include { VEP } from '../modules/vep'
include { ANNOTATE } from '../modules/annotate'

workflow VARIANTS {
    
    take:
    chrm_bam
    hmtdb
    numts
    fasta
    vep_cache
    gnomad_vcf

    main:

    // if the data is rna-seq, skip numt filtering
    if ( "${params.rnaseq}" == "false" ) {
        
        RTN (
            chrm_bam,
            hmtdb,
            numts
            )
        
        pileup_input = RTN.out.numt_filt_ch
        
    } else {
        
        pileup_input = chrm_bam
            .map({ meta, bam, index -> [meta, bam] })

    }

    PILEUP (
        pileup_input,
        fasta
    )

    if ( "${params.matched}" == "true" ) {
        
        // Split pileup files in tumour and normal samples
        ch_tumour = PILEUP.out.mpileup_ch
            .filter({ meta, pileup -> meta.group == 'tumour' })
            .map({ meta, pileup -> [[sample: meta.sample], pileup] })

        ch_normal = PILEUP.out.mpileup_ch
            .filter({ meta, pileup -> meta.group == 'normal' })
            .map({ meta, pileup -> [[sample: meta.sample], pileup] })
        
        // Join pileup files together based on sample ID 
        matched_pileup_ch = ch_tumour
            .join(ch_normal)

    } else {

        matched_pileup_ch = PILEUP.out.mpileup_ch
            .map({ meta, pileup -> [[sample: meta.sample], pileup, file("$baseDir/misc/main.nf")] })

    }

    VARSCAN2 (
        matched_pileup_ch
    )

    COMBINE_VARIANTS (
        VARSCAN2.out.indels_ch.collect(),
        VARSCAN2.out.snvs_ch.collect()
    )

    VEP_FORMAT (
        COMBINE_VARIANTS.out.combined_ch
    )

    VEP (
        VEP_FORMAT.out.vep_formatted_ch,
        fasta,
        vep_cache
    )

    ANNOTATE (
        COMBINE_VARIANTS.out.combined_ch,
        VEP.out.annotated_ch,
        gnomad_vcf
    )

}