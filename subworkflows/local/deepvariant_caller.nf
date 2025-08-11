//
// Call variants with Deepvariant
//

include { DEEPVARIANT_RUNDEEPVARIANT as DEEPVARIANT       }   from '../../modules/nf-core/deepvariant/rundeepvariant/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF          }   from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_GVCF         }   from '../../modules/nf-core/bcftools/concat/main'
include { DEEPVARIANT_VCFSTATSREPORT as VCF_STATS_REPORT  }   from '../../modules/nf-core/deepvariant/vcfstatsreport/main'
include { DEEPVARIANT_VCFSTATSREPORT as GVCF_STATS_REPORT }   from '../../modules/nf-core/deepvariant/vcfstatsreport/main'
include { TABIX_BGZIP as BGZIP                            }   from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as TABIX_CSI                        }   from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TBI                        }   from '../../modules/nf-core/tabix/tabix/main'

workflow DEEPVARIANT_CALLER {
    take:
    reads_fasta    // [ val(meta), cram, crai, interval, fasta_file_name, fasta, fai ]
    max_length     // [ val(max_length) - maximum chromosome length in the fasta file  ]

    main:
    ch_versions = Channel.empty()

    reads_fasta.map { meta, cram, crai, interval, fasta_file_name, fasta, fai ->
                     [ [ id: meta.id + "_" + fasta_file_name,
                         sample: meta.id,
                         type: meta.datatype,
                         fasta_file_name: fasta_file_name
                       ],
                       cram,
                       crai,
                       interval
                     ] }
               .set { cram_crai }

    // fasta
    fasta = reads_fasta.map { meta, cram, crai, interval, fasta_file_name, fasta, fai ->
                             [ [
                                id: meta.id + "_" + fasta_file_name,
                                sample: meta.id,
                                type: meta.datatype ],
                              fasta
                             ]
                            }

    // fai
    fai = reads_fasta.map{ meta, cram, crai, interval, fasta_file_name, fasta, fai ->
                           [ [ id: meta.id + "_" + fasta_file_name, sample: meta.id, type: meta.datatype ],
                             fai
                           ]
                         }

    // split fasta in compressed format, no gzi index file needed
    gzi = [ [], [] ]
    par_bed = [ [], [] ]

    // call deepvariant
    DEEPVARIANT ( cram_crai, fasta, fai, gzi, par_bed )
    ch_versions = ch_versions.mix ( DEEPVARIANT.out.versions.first() )

    // group the vcf files together by sample
    DEEPVARIANT.out.vcf
        .join(DEEPVARIANT.out.vcf_index)
        .map { meta, vcf, index -> [
            [ id: meta.fasta_file_name.tokenize(".")[0..-2].join(".")
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
            [ id: meta.fasta_file_name.tokenize(".")[0..-2].join(".")
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

    // compress the vcf and gvcf files
    vcf_to_compress   = BCFTOOLS_CONCAT_VCF.out.vcf.mix ( BCFTOOLS_CONCAT_GVCF.out.vcf )
    ch_compressed_vcf = BGZIP ( vcf_to_compress ).output
    ch_versions       = ch_versions.mix ( BGZIP.out.versions.first() )

    // index the compressed files in two formats for maximum compatibility (each has its own limitation)
    // select the type of index to use based on the maximum sequence length
    ch_compressed_vcf
        .combine(max_length)
        .map { meta_vcf, vcf, meta -> [ meta_vcf + meta, vcf ] }
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

    // generate vcf stats report
    VCF_STATS_REPORT ( BCFTOOLS_CONCAT_VCF.out.vcf )
    ch_versions = ch_versions.mix ( VCF_STATS_REPORT.out.versions.first() )

    // generate g vcf stats report
    GVCF_STATS_REPORT ( BCFTOOLS_CONCAT_GVCF.out.vcf )
    ch_versions = ch_versions.mix ( GVCF_STATS_REPORT.out.versions.first() )

    emit:
    vcf      = BCFTOOLS_CONCAT_VCF.out.vcf           // channel: [ val(meta), path(vcf) ]
    gvcf     = BCFTOOLS_CONCAT_GVCF.out.vcf          // channel: [ val(meta), path(gvcf) ]
    compressed_vcf    = ch_compressed_vcf            // channel: [ val(meta), path(output)]
    vcf_csi  = ch_indexed_vcf_csi                    // channel: [ val(meta), path(csi)]
    vcf_tbi  = ch_indexed_vcf_tbi                    // channel: [ val(meta), path(tbi)]
    vcf_stats_report  = VCF_STATS_REPORT.out.report  // channel: [ val(meta), path(report) ]
    gvcf_stats_report = GVCF_STATS_REPORT.out.report // channel: [ val(meta), path(report) ]
    versions = ch_versions                           // channel: [ versions.yml ]
}
