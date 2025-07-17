//
// Run Himut to analyse VCF
//

include { HIMUT                    }   from '../../modules/local/himut.nf'
include { TABIX_BGZIP as BGZIP     }   from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as TABIX_CSI }   from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TBI }   from '../../modules/nf-core/tabix/tabix/main'

workflow RUN_HIMUT {
    take:
    fasta                // [ val(meta), fasta           ]
    fasta_index          // [ val(meta), fai             ]
    assembly_report      // [ val(meta), assembly_report ]
    bam                  // [ val(meta), bam             ]
    bam_index            // [ val(meta), bai             ]
    vcf_input            // [ val(meta), vcf_input       ]
    vcf_index            // [ val(meta), vcf_tbi         ]
    max_length           // [ val(max_length)            ]

    main:
    ch_versions = Channel.empty()

    // run Himut
    HIMUT ( fasta, fasta_index, assembly_report, bam, bam_index, vcf_input, vcf_index )
    ch_versions = ch_versions.mix ( HIMUT.out.versions.first() )

    // compress the vcf outputs of Himut
    himut_vcf_to_compress   = HIMUT.out.vcf_output.mix ( HIMUT.out.smm_vcf_output )
    ch_compressed_himut_vcf = BGZIP ( himut_vcf_to_compress ).output
    ch_versions             = ch_versions.mix ( BGZIP.out.versions.first() )


    // index the compressed himut vcf files in two formats for maximum compatibility (each has its own limitation)
    // select the type of index to use based on the maximum sequence length
    ch_compressed_himut_vcf
        .combine( max_length )
        .map { meta_vcf, vcf, meta -> [ meta_vcf + meta, vcf ] }
        .branch { meta, vcf ->
        tbi_and_csi: meta.max_length < 2**29
        only_csi:    meta.max_length < 2**32
        }
        .set { himut_tabix_selector }

    // do the indexing on the compatible gvcf files
    ch_himut_vcf_csi   = TABIX_CSI ( himut_tabix_selector.tbi_and_csi.mix(himut_tabix_selector.only_csi) ).csi
    ch_versions        = ch_versions.mix ( TABIX_CSI.out.versions.first() )
    ch_himut_vcf_tbi   = TABIX_TBI ( himut_tabix_selector.tbi_and_csi ).tbi
    ch_versions        = ch_versions.mix ( TABIX_TBI.out.versions.first() )

    emit:
    himut_vcf            = HIMUT.out.vcf_output     // channel: [ val(meta), path(vcf_output)     ]
    himut_smm_vcf        = HIMUT.out.smm_vcf_output // channel: [ val(meta), path(smm_vcf_output) ]
    compressed_himut_vcf = ch_compressed_himut_vcf  // channel: [ val(meta), path(output)         ]
    himut_vcf_csi        = ch_himut_vcf_csi         // channel: [ val(meta), path(csi)            ]
    himut_vcf_tbi        = ch_himut_vcf_tbi         // channel: [ val(meta), path(tbi)            ]
    versions = ch_versions
}
