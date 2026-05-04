//
// Call variants with Deepvariant
//

include { DEEPVARIANT_RUNDEEPVARIANT as DEEPVARIANT       }   from '../../modules/nf-core/deepvariant/rundeepvariant/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF          }   from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_GVCF         }   from '../../modules/nf-core/bcftools/concat/main'
include { TABIX_BGZIP as BGZIP                            }   from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as TABIX_CSI                        }   from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TBI                        }   from '../../modules/nf-core/tabix/tabix/main'

workflow DEEPVARIANT_CALLER {
    take:
    reads_fasta    // [ val(meta), bam, bai, interval, val(meta_fasta), fasta, fai ]
    max_length     // [ val(meta_max_length) - maximum chromosome length in the fasta file  ]

    main:
    ch_versions = Channel.empty()

    reads_fasta.map { meta, bam, bai, interval, meta_fasta, fasta, fai ->
                     [ [ id: meta.id + "_" + meta_fasta.id,
                         sample: meta.id,
                         type: meta.datatype,
                         fasta_id: meta_fasta.id.tokenize(".")[0..-2].join(".") // Strip the suffix added by seqkit
                       ],
                       bam,
                       bai,
                       interval
                     ] }
               .set { bam_bai }

    // fasta
    fasta = reads_fasta.map { meta, bam, bai, interval, meta_fasta, fasta, fai -> [ meta_fasta, fasta ] }

    // fai
    fai = reads_fasta.map{ meta, bam, bai, interval, meta_fasta, fasta, fai -> [ meta_fasta, fai ] }

    // split fasta in compressed format, no gzi index file needed
    gzi = [ [], [] ]
    par_bed = [ [], [] ]

    // call deepvariant
    DEEPVARIANT ( bam_bai, fasta, fai, gzi, par_bed )
    ch_versions = ch_versions.mix ( DEEPVARIANT.out.versions.first() )

    // group the vcf files together by sample
    DEEPVARIANT.out.vcf
        .join(DEEPVARIANT.out.vcf_index)
        .map { meta, vcf, index -> [
            [ id: meta.fasta_id
                + "." + meta.type
                + "." + meta.sample
            ],
            vcf,
            index
        ] }
        .groupTuple()
        .set { vcf }

    // concat vcf files
    BCFTOOLS_CONCAT_VCF ( vcf )
    ch_versions = ch_versions.mix ( BCFTOOLS_CONCAT_VCF.out.versions.first() )

    // group the g vcf files together by sample
    DEEPVARIANT.out.gvcf
        .join(DEEPVARIANT.out.gvcf_index)
        .map { meta, gvcf, index -> [
            [ id: meta.fasta_id
                + "." + meta.type
                + "." + meta.sample
            ],
            gvcf,
            index
        ] }
        .groupTuple()
        .set { g_vcf }

    // concat g vcf files
    BCFTOOLS_CONCAT_GVCF ( g_vcf )
    ch_versions = ch_versions.mix ( BCFTOOLS_CONCAT_GVCF.out.versions.first() )

    // we'll want to index the vcf in two formats for maximum compatibility (each has its own limitation)
    // selection of the type of index is based on the maximum sequence length
    BCFTOOLS_CONCAT_VCF.out.vcf
        .combine( max_length )
        .map { meta_vcf, vcf, meta -> [ meta + meta_vcf, vcf ] }
        .set { ch_vcf_with_seq_lengths }

    // compress the vcf file as the input of Himut
    ch_compressed_vcf = BGZIP ( ch_vcf_with_seq_lengths ).output
    ch_versions       = ch_versions.mix ( BGZIP.out.versions.first() )

    // do the selection
    ch_compressed_vcf
        .branch { meta, vcf ->
        tbi_and_csi: meta.max_length < 2**29
        only_csi:    meta.max_length < 2**32
        }
        .set { tabix_selector }

    // do the indexing on the compatible gvcf files
    ch_indexed_vcf_csi = TABIX_CSI ( tabix_selector.tbi_and_csi.mix(tabix_selector.only_csi) ).csi
    ch_versions        = ch_versions.mix ( TABIX_CSI.out.versions.first() )
    ch_indexed_vcf_tbi = TABIX_TBI ( tabix_selector.tbi_and_csi ).tbi
    ch_versions        = ch_versions.mix ( TABIX_TBI.out.versions.first() )

    emit:
    vcf      = BCFTOOLS_CONCAT_VCF.out.vcf           // channel: [ val(meta), path(vcf)    ]
    gvcf     = BCFTOOLS_CONCAT_GVCF.out.vcf          // channel: [ val(meta), path(gvcf)   ]
    compressed_vcf    = ch_compressed_vcf            // channel: [ val(meta), path(output) ]
    vcf_csi  = ch_indexed_vcf_csi                    // channel: [ val(meta), path(csi)    ]
    vcf_tbi  = ch_indexed_vcf_tbi                    // channel: [ val(meta), path(tbi)    ]
    versions = ch_versions                           // channel: [ versions.yml            ]
}
