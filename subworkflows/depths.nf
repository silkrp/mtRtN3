//
// Get sequencing depths for each base-pair on the mtDNA
//

include { MTDNA_DEPTHS } from '../modules/mtdna_depths'
include { COMBINE_DEPTHS } from '../modules/combine_depths'

workflow DEPTHS {
    
    take:
    chrm_bam

    main:

    MTDNA_DEPTHS (
        chrm_bam
    )

    COMBINE_DEPTHS (
        MTDNA_DEPTHS.out.depths_ch.collect()
    )

}