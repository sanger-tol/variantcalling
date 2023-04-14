//
// Split input fasta file by sequence and filter the input reads
//

include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main'

workflow INPUT_FILTER_SPLIT {

    take:
    fasta    // path to reference fasta file, either compressed or uncompressed
    fai      // path to compressed or uncompressed reference fasta index file
    gzi      // path to compressed fasta index file

    reads    // [ val(meta), data, index ]
    interval // path to interval bed file

    main:
    ch_versions = Channel.empty()

    // split the fasta file into files with one sequence each
    Channel
     .fromPath(fasta)
     .splitFasta(file:true)
     .map{ [ [id: 'splitting_fasta'],  it ] }
     .set{ splitted_fasta }

    // index splitted fasta files
    SAMTOOLS_FAIDX ( splitted_fasta  )
    ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions.first() )


    // filter reads
    SAMTOOLS_VIEW ( reads, fasta, [] )
    ch_versions = ch_versions.mix ( SAMTOOLS_VIEW.out.versions.first() )
    
    // combine reads with splitted references
    cram_crai_fasta_fai = SAMTOOLS_VIEW.out.cram
                        .join(SAMTOOLS_VIEW.out.crai)
                        .map { filtered_reads -> filtered_reads + [interval ?: []] }
                        .combine( splitted_fasta.join( SAMTOOLS_FAIDX.out.fai ) )


    emit:
    reads_fasta    = cram_crai_fasta_fai    // channel:[ val(meta), cram, crai, interval, val(meta), fasta, fai ]
    versions       = ch_versions            // channel: [ versions.yml ]

}
