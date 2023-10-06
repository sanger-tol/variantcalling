/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowVariantcalling.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.fasta, params.fai, params.gzi, params.interval ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = Channel.fromPath(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = Channel.fromPath(params.fasta) } else { exit 1, 'Reference fasta not specified!'   }

// Check optional parameters
if (params.fai)     { ch_fai   = Channel.fromPath(params.fai)         } else { ch_fai      = Channel.empty() }
if (params.gzi)     { ch_gzi   = Channel.fromPath(params.gzi)         } else { ch_gzi      = Channel.empty() }
if (params.interval){ ch_interval = Channel.fromPath(params.interval) } else { ch_interval = Channel.empty() }

if (params.sort_input)          { sort_input    = params.sort_input              } else { sort_input    = false       }
if (params.split_fasta_cutoff ) { split_fasta_cutoff = params.split_fasta_cutoff } else { split_fasta_cutoff = 100000 }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK        } from '../subworkflows/local/input_check'
include { INPUT_MERGE        } from '../subworkflows/local/input_merge'
include { INPUT_FILTER_SPLIT } from '../subworkflows/local/input_filter_split'
include { DEEPVARIANT_CALLER } from '../subworkflows/local/deepvariant_caller'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow VARIANTCALLING {

    ch_versions = Channel.empty()

    //
    // check reference fasta index given or not
    //
    if( params.fai == null || ( params.fasta.endsWith('fasta.gz') && params.gzi == null ) ){ 
   
       ch_fasta
        .map { fasta -> [ [ 'id': fasta.baseName ], fasta ] }
        .set { ch_genome }

       SAMTOOLS_FAIDX ( ch_genome,  [[], []])
       ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )
       
       SAMTOOLS_FAIDX.out.fai
        .map{ mata, fai -> fai }
        .set{ ch_fai }
  
       SAMTOOLS_FAIDX.out.gzi
        .map{ meta, gzi -> gzi }
        .set{ ch_gzi }

    }

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: merge the input reads by sample name
    //
    INPUT_MERGE (
        ch_fasta,
        ch_fai,
        ch_gzi,
        INPUT_CHECK.out.reads,
        sort_input
    )
    ch_versions = ch_versions.mix(INPUT_MERGE.out.versions)


    //
    // SUBWORKFLOW: split the input fasta file and filter input reads
    //
    INPUT_FILTER_SPLIT (
        ch_fasta,
        ch_fai,
        ch_gzi,
        INPUT_MERGE.out.indexed_merged_reads,
        ch_interval,
        split_fasta_cutoff
    )
    ch_versions = ch_versions.mix(INPUT_FILTER_SPLIT.out.versions)

    //
    // SUBWORKFLOW: call deepvariant
    //
    DEEPVARIANT_CALLER (
        INPUT_FILTER_SPLIT.out.reads_fasta
    )
    ch_versions = ch_versions.mix(DEEPVARIANT_CALLER.out.versions)

    //
    // MODULE: Combine different version together
    // 
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

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
