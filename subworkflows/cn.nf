//
// Get sequencing depths for each base-pair on the mtDNA
//

include { COVERAGE } from '../modules/coverage'
include { CALCULATE_CN } from '../modules/calculate_cn'
include { COMBINE_CN } from '../modules/combine_cn'

workflow CN {
    
    take:
    bams
    cram_fasta

    main:

    COVERAGE (
        bams,
        cram_fasta
    )

    CALCULATE_CN (
        COVERAGE.out.coverage_ch
    )

    COMBINE_CN (
        CALCULATE_CN.out.cn_ch.collect()
    )

}