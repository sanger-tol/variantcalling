/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALIGN_PACBIO       } from '../subworkflows/local/align_pacbio'
include { INPUT_MERGE        } from '../subworkflows/local/input_merge'
include { INPUT_FILTER_SPLIT } from '../subworkflows/local/input_filter_split'
include { DEEPVARIANT_CALLER } from '../subworkflows/local/deepvariant_caller'
include { PROCESS_VCF        } from '../subworkflows/local/process_vcf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_variantcalling_pipeline'
include { SAMTOOLS_FAIDX         } from '../modules/nf-core/samtools/faidx/main'
include { UNTAR                  } from '../modules/nf-core/untar/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow VARIANTCALLING {

    take:
    ch_reads            // channel: samplesheet read in from --input
    ch_fasta            // channel: fasta file read in from --fasta
    ch_fai              // channel: fai file read in from --fai
    ch_interval         // channel: interval file read in from --interval
    ch_positions        // channel: positions to include or exclude in the variant calling

    main:
    ch_versions = Channel.empty()

    //
    // channel for reference genome
    //
    // Remenber to fix the fasta.size with total_length in the next merge
    // Also remember to change from toInterger to toLong
    ch_fasta
        .map { fasta -> [ [ 'id': fasta.baseName -  ~/.fa\w*$/ , 'genome_size': fasta.size() ], fasta ] }
        .first()
        .set { ch_genome }

    // reference genome parameter
    if (params.split_fasta_cutoff ) { split_fasta_cutoff = params.split_fasta_cutoff } else { split_fasta_cutoff = 100000 }

    //
    // check reference fasta index given or not
    //
    if( !params.fai ) {

        SAMTOOLS_FAIDX ( ch_genome,  [[], []] )
        ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )

        // generate fai that is used to determine the maximum length of chromosome
        ch_genome_index_fai = SAMTOOLS_FAIDX.out.fai
        ch_genome_index     = params.fasta.endsWith('.gz') ? SAMTOOLS_FAIDX.out.gzi : SAMTOOLS_FAIDX.out.fai

    } else {
        ch_fai
            .map { fai -> [ [ 'id': fai.baseName ], fai ] }
            .first()
            .set { ch_genome_index }

        ch_genome_index_fai  = ch_genome_index
        if ( !params.fai.endsWith(".fai") ) {
            ch_genome_index_fai = SAMTOOLS_FAIDX ( ch_genome,  [[], []] ).fai
            ch_versions         = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )
        }
    }

    ch_genome_index_fai
        .map { meta, fai_file -> [ [ id: meta.id ] + get_sequence_map(fai_file) ] }
        .set { ch_genome_info }

    //
    // SUBWORKFLOW: align reads if required
    //
    if( params.align ) {

        if ( params.vector_db.endsWith( '.tar.gz' ) ) {

            UNTAR ( [ [:], params.vector_db ] ).untar
            | map { meta, file -> file }
            | set { ch_vector_db }
            ch_versions = ch_versions.mix ( UNTAR.out.versions )


        } else {

            Channel.fromPath ( params.vector_db )
            | set { ch_vector_db }

        }

        ALIGN_PACBIO (
            ch_genome,
            ch_reads,
            ch_vector_db
        )
        ch_versions = ch_versions.mix( ALIGN_PACBIO.out.versions )

        ALIGN_PACBIO.out.cram
            .join( ALIGN_PACBIO.out.crai )
            .set{ ch_aligned_reads }

    } else {

        //
        // SUBWORKFLOW: merge the input reads by sample name
        //
        INPUT_MERGE (
            ch_genome,
            ch_genome_index,
            ch_reads
        )
        ch_versions = ch_versions.mix( INPUT_MERGE.out.versions )
        ch_aligned_reads = INPUT_MERGE.out.indexed_merged_reads

    }

    //
    // SUBWORKFLOW: split the input fasta file and filter input reads
    //
    INPUT_FILTER_SPLIT (
        ch_fasta,
        ch_aligned_reads,
        ch_interval,
        split_fasta_cutoff
    )
    ch_versions = ch_versions.mix( INPUT_FILTER_SPLIT.out.versions )


    //
    // SUBWORKFLOW: call deepvariant
    //
    DEEPVARIANT_CALLER (
        INPUT_FILTER_SPLIT.out.reads_fasta,
        ch_genome_info
    )
    ch_versions = ch_versions.mix( DEEPVARIANT_CALLER.out.versions )


    //
    // convert VCF channel meta id
    //
    DEEPVARIANT_CALLER.out.vcf
        .map{ meta, vcf -> [ [ id: vcf.baseName ], vcf ] }
        .set{ vcf }

    //
    // process VCF output files
    //
    PROCESS_VCF( vcf, ch_positions )
    ch_versions = ch_versions.mix( PROCESS_VCF.out.versions )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'variantcalling_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}


//
// FUNCTION: get sequence map
// Read the .fai file, extract sequence statistics, and make an extended meta map
//

def get_sequence_map(fai_file) {
    def n_sequences    = 0
    def max_length     = 0
    def total_length   = 0
    fai_file.eachLine { line ->
        def lspl       = line.split('\t')
        def chrom      = lspl[0]
        def length     = lspl[1].toInteger()
        n_sequences ++
        total_length  += length
        if (length > max_length) {
            max_length = length
        }
    }

    def sequence_map = [:]
    sequence_map.n_sequences    = n_sequences
    sequence_map.total_length   = total_length
    if (n_sequences) {
        sequence_map.max_length = max_length
    }
    return sequence_map
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
