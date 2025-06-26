//
// Check input samplesheet and get read channels
//

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    Channel
        .fromPath( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { create_bam_channel(it) }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]

}

// Function to get list of [ meta, [ bam, index ] ]
def create_bam_channel(LinkedHashMap row) {
    
    // create meta map
    def meta = [:]
    meta.sample = row.sample
    meta.group  = row.group

    // check if bam and index files exist
    def bam_meta = []
    if (!file(row.bam).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> BAM file does not exist!\n${row.bam}"
    }
    if (!file(row.index).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> index file does not exist!\n${row.bam}"
    }

    // add path(s) of the bam file(s) to the meta map
    bam_meta = [ meta, file(row.bam), file(row.index) ]
    return bam_meta

}