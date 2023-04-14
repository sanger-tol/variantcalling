//
// call deepvariant
//

include { DEEPVARIANT }   from '../../modules/nf-core/deepvariant/main'

workflow DEEPVARIANT_CALLER {

    take:
    reads_fasta    // [ val(meta), cram, crai, interval, val(meta), fasta, fai ]

    main:
    ch_versions = Channel.empty()

    //DEEPVARIANT( cram_crai, fasta, fai, gzi )
    //ch_versions = ch_versions.mix ( DEEPVARIANT.out.versions )


    emit:
    //vcf      = DEEPVARIANT.out.vcf               // /path/to/vcf
    //gvcf      = DEEPVARIANT.out.gvcf             // /path/to/gvcf
    versions = ch_versions                       // channel: [ versions.yml ]
}
