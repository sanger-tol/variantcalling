//
// Check input samplesheet and get read channels
//

include { SAMTOOLS_VIEW } from '../../modules/nf-core/samtools/view/main'
include { DEEPVARIANT }   from '../../modules/nf-core/deepvariant/main'

workflow DEEPVARIANT_CALLER {

    take:
    reads    // [ val(meta), data, index ]
    fasta    // path to reference fasta file, either compressed or uncompressed
    fai      // path to compressed or uncompressed reference fasta index file
    gzi      // path to compressed fasta index file
    interval // path to interval bed file

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_VIEW ( reads, fasta, [] )
    ch_versions = ch_versions.mix ( SAMTOOLS_VIEW.out.versions )
     
    cram_crai = SAMTOOLS_VIEW.out.cram
                .join(SAMTOOLS_VIEW.out.crai)
                .map { it + [interval ?: []] }

    if ( !gzi ) {
        gzi = []
    }

    DEEPVARIANT( cram_crai, fasta, fai, gzi )
    ch_versions = ch_versions.mix ( DEEPVARIANT.out.versions )


    emit:
    cram     = SAMTOOLS_VIEW.out.cram            // channel: [ val(meta), /path/to/cram ]
    crai     = SAMTOOLS_VIEW.out.crai            // channel: [ val(meta), /path/to/crai ]
    vcf      = DEEPVARIANT.out.vcf               // /path/to/vcf
    gvcf      = DEEPVARIANT.out.gvcf             // /path/to/gvcf
    versions = ch_versions                       // channel: [ versions.yml ]
}
