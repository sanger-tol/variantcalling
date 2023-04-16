//
// call deepvariant
//

include { DEEPVARIANT }   from '../../modules/nf-core/deepvariant/main'

workflow DEEPVARIANT_CALLER {

    take:
    reads_fasta    // [ val(meta), cram, crai, interval, fasta_file_name, fasta, fai ]

    main:
    ch_versions = Channel.empty()

    // [ val(meta), cram, crai, interval ]
    reads_fasta.map{ 
                    [ [ id: it[0].id + "_" + it[4], sample: it[0].id, type: it[0].type ], 
                       it[1],
                       it[2], 
                       it[3] 
                    ] }
               .set{ cram_crai }

    // fasta
    reads_fasta.map{ [ it[5] ] }
            .set{ fasta }

    // fai
    reads_fasta.map{ [ it[6] ] }
            .set{ fai }

    // split fasta in compressed format, no gzi index file needed
    gzi = []

    // call deepvariant
    DEEPVARIANT( cram_crai, fasta, fai, gzi )
    ch_versions = ch_versions.mix ( DEEPVARIANT.out.versions.first() )

    emit:
    vcf      = DEEPVARIANT.out.vcf              // /path/to/vcf
    gvcf     = DEEPVARIANT.out.gvcf             // /path/to/gvcf
    versions = ch_versions                      // channel: [ versions.yml ]
}
