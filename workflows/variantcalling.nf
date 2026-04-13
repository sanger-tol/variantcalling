/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALIGN_PACBIO           } from '../subworkflows/local/align_pacbio'
include { INPUT_MERGE            } from '../subworkflows/local/input_merge'
include { INPUT_FILTER_SPLIT     } from '../subworkflows/local/input_filter_split'
include { DEEPVARIANT_CALLER     } from '../subworkflows/local/deepvariant_caller'

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
    ch_reads // channel: samplesheet read in from --input
    ch_fasta // channel: fasta file read in from --fasta
    ch_interval // channel: interval file read in from --interval

    main:
    ch_versions = channel.empty()

    //
    // Channel for reference genome
    //
    // Remenber to fix the fasta.size with total_length in the next merge
    ch_genome = ch_fasta.map { fasta ->
        [
            [
                'id': fasta.baseName - ~/.fa\w*$/,
                'single_end': true,
            ],
            fasta,
        ]
    }


    SAMTOOLS_FAIDX(ch_genome, [[], []])
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // generate fai that is used to determine the maximum length of chromosome
    // also add the gzi if present as it is needed for bgzip-ed genomes
    ch_genome_info = ch_genome
       .join( SAMTOOLS_FAIDX.out.fai )
       .join( SAMTOOLS_FAIDX.out.gzi, remainder: true )
       .map { meta, fa, fai, gzi ->
           def index_file = (fa.name.endsWith('.gz') && gzi) ? [fai, gzi] : fai
           [meta + get_sequence_map(fai), fa, index_file]
        }
        .collect()
        .multiMap { meta, fa, fai ->
            meta: meta
            fasta: [meta, fa]
            index: [meta, fai]
        }


    //
    // SUBWORKFLOW: align reads if required
    //
    if (params.align) {

        if (params.vector_db.endsWith('.tar.gz')) {

            ch_vector_db = UNTAR([[:], params.vector_db]).untar.map { _meta, file -> file }
            ch_versions = ch_versions.mix(UNTAR.out.versions)
        }
        else {

            ch_vector_db = channel.fromPath(params.vector_db)
        }

        ALIGN_PACBIO(
            ch_genome_info.fasta,
            ch_reads,
            ch_vector_db,
        )
        ch_versions = ch_versions.mix(ALIGN_PACBIO.out.versions)

        ch_aligned_reads = ALIGN_PACBIO.out.cram.join(ALIGN_PACBIO.out.crai)
    }
    else {

        //
        // SUBWORKFLOW: merge the input reads by sample name
        //
        INPUT_MERGE(
            ch_genome_info.fasta,
            ch_genome_info.index,
            ch_reads,
        )
        ch_versions = ch_versions.mix(INPUT_MERGE.out.versions)
        ch_aligned_reads = INPUT_MERGE.out.indexed_merged_reads
    }


    //
    // SUBWORKFLOW: split the input fasta file and filter input reads
    //
    INPUT_FILTER_SPLIT(
        ch_genome_info.fasta,
        ch_aligned_reads,
        ch_interval,
    )
    ch_versions = ch_versions.mix(INPUT_FILTER_SPLIT.out.versions)


    //
    // SUBWORKFLOW: call deepvariant
    //
    DEEPVARIANT_CALLER(
        INPUT_FILTER_SPLIT.out.reads_fasta,
        ch_genome_info.meta.max_length,
    )


    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'variantcalling_software_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    emit:
    versions = ch_collated_versions // channel: [ path(versions.yml) ]
}


//
// FUNCTION: get sequence map
// Read the .fai file, extract sequence statistics, and make an extended meta map
//

def get_sequence_map(fai_file) {
    def n_sequences = 0
    def max_length = 0
    def total_length = 0
    fai_file.eachLine { line ->
        def lspl = line.split('\t')
        // def chrom      = lspl[0]
        def length = lspl[1].toLong()
        n_sequences += 1
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
