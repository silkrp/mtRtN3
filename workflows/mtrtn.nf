/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// NEED TO HAVE SOME CHECK TO MAKE SURE EVERYTHING IS SET TO A SENSIBLE VALUE E.G. TRUE / FALSE

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES
//
include { PARSE_MTDNA } from '../modules/parse_mtdna'

//
// SUBWORKFLOWS
//
include { PREPARE_GENOME } from '../subworkflows/prepare_genome'
include { INPUT_CHECK } from '../subworkflows/input_check'
include { VARIANTS } from '../subworkflows/variants'
include { DEPTHS } from '../subworkflows/depths'
include { CN } from '../subworkflows/cn'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MTRTN {

    //
    // SUBWORKFLOW: Prepare reference genome files
    // 
    PREPARE_GENOME (
        params.fasta,
        params.hmtdb,
        params.numts,
        params.cram_fasta,
        params.vep_cache,
        params.vep_species,
        params.vep_genome,
        params.vep_cache_version,
        params.gnomad_vcf
    ) 

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .set { ch_bam } 

    //
    // MODULE: Parse mtDNA reads from alignment file using samtools 
    //
    PARSE_MTDNA (
        ch_bam,
        PREPARE_GENOME.out.cram_fasta
    )

    //
    // Parse subworkflows and run each specified
    //
    def selected_workflows = params.workflows ? params.workflows.split(',') : []

    if ('variant-calling' in selected_workflows) {
        
        //
        // SUBWORKFLOW: mtDNA variant calling
        // 
        VARIANTS (
            PARSE_MTDNA.out.chrm_bam,
            PREPARE_GENOME.out.hmtdb,
            PREPARE_GENOME.out.numts,
            PREPARE_GENOME.out.fasta,
            PREPARE_GENOME.out.vep_cache,
            PREPARE_GENOME.out.gnomad_vcf
        )
    }

    if ('cn' in selected_workflows) {
        
        //
        // SUBWORKFLOW: mtDNA copy-number
        // 
        CN (
            ch_bam,
            PREPARE_GENOME.out.cram_fasta
        )
    }

    if ('depths' in selected_workflows) {
        
        //
        // SUBWORKFLOW: mtDNA sequencing depths
        // 
        DEPTHS (
            PARSE_MTDNA.out.chrm_bam
        )
    }

}

