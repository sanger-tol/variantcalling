/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowVariantcalling.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.fasta, params.fai, params.interval, params.include_positions, params.exclude_positions ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = Channel.fromPath(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = Channel.fromPath(params.fasta) } else { exit 1, 'Reference fasta not specified!'   }

// Check optional parameters
if (params.fai){
    if( ( params.fasta.endsWith('.gz') && params.fai.endsWith('.fai') )
        ||
        ( !params.fasta.endsWith('.gz') && params.fai.endsWith('.gzi') )
    ){
        exit 1, 'Reference fasta and its index file format not matched!'
    }
    ch_fai = Channel.fromPath(params.fai)
} else {
    ch_fai = Channel.empty()
}

if (params.interval){ ch_interval = Channel.fromPath(params.interval) } else { ch_interval = Channel.empty() }

if (params.split_fasta_cutoff ) { split_fasta_cutoff = params.split_fasta_cutoff } else { split_fasta_cutoff = 100000 }

if ( (params.include_positions) && (params.exclude_positions) ){
    exit 1, 'Only one positions file can be given to include or exclude!'
}else if (params.include_positions){
    ch_positions = Channel.fromPath(params.include_positions)
} else if (params.exclude_positions){
    ch_positions = Channel.fromPath(params.exclude_positions)
} else {
    ch_positions = []
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK        } from '../subworkflows/local/input_check'
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

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SAMTOOLS_FAIDX              } from '../modules/nf-core/samtools/faidx/main'
include { UNTAR                       } from '../modules/nf-core/untar/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow VARIANTCALLING {

    ch_versions = Channel.empty()
    ch_fasta
        .map { fasta -> [ [ 'id': fasta.baseName -  ~/.fa\w*$/ , 'genome_size': fasta.size() ], fasta ] }
        .first()
        .set { ch_genome }

    //
    // check reference fasta index given or not
    //

    if( params.fai == null ){

        SAMTOOLS_FAIDX ( ch_genome,  [[], []] )
        ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )

        // generate fai that is used to determine the maximum length of chromosome
        ch_genome_index_fai = SAMTOOLS_FAIDX.out.fai
        ch_genome_index = params.fasta.endsWith('.gz') ? SAMTOOLS_FAIDX.out.gzi : SAMTOOLS_FAIDX.out.fai

    }else{
        ch_fai
            .map { fai -> [ [ 'id': fai.baseName ], fai ] }
            .first()
            .set { ch_genome_index }

        ch_genome_index_fai  = ch_genome_index
        if ( !params.fai.endsWith(".fai") ) {
            ch_genome_index_fai = SAMTOOLS_FAIDX ( ch_genome,  [[], []] ).fai
            ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )
        }
    }

    ch_genome_index_fai
        .map { meta, index -> [ [ id: meta.id ] + get_sequence_map(index) ] }
        .set { ch_genome_info }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix( INPUT_CHECK.out.versions )


    //
    // SUBWORKFLOW: align reads if required
    //
    if( params.align ){

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
            INPUT_CHECK.out.reads,
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
            INPUT_CHECK.out.reads,
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
    // MODULE: Combine different version together
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}


//
// function: get sequence map
//

// Read the .fai file, extract sequence statistics, and make an extended meta map
def get_sequence_map(fai_file) {
    def n_sequences = 0
    def max_length = 0
    def total_length = 0
    fai_file.eachLine { line ->
        def lspl   = line.split('\t')
        def chrom  = lspl[0]
        def length = lspl[1].toInteger()
        n_sequences ++
        total_length += length
        if (length > max_length) {
            max_length = length
        }
    }

    def sequence_map = [:]
    sequence_map.n_sequences = n_sequences
    sequence_map.total_length = total_length
    if (n_sequences) {
        sequence_map.max_length = max_length
    }
    return sequence_map
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
