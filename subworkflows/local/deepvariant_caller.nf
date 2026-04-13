//
// Call variants with Deepvariant
//

include { DEEPVARIANT_RUNDEEPVARIANT as DEEPVARIANT       } from '../../modules/nf-core/deepvariant/rundeepvariant/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF          } from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_GVCF         } from '../../modules/nf-core/bcftools/concat/main'
include { DEEPVARIANT_VCFSTATSREPORT as VCF_STATS_REPORT  } from '../../modules/nf-core/deepvariant/vcfstatsreport/main'
include { DEEPVARIANT_VCFSTATSREPORT as GVCF_STATS_REPORT } from '../../modules/nf-core/deepvariant/vcfstatsreport/main'
include { BGZIPTABIX                                      } from '../../modules/sanger-tol/bgziptabix/main'

workflow DEEPVARIANT_CALLER {
    take:
    reads_fasta // [ val(meta), cram, crai, interval, val(meta_fasta), fasta, fai ]
    max_length // [ val(meta_max_length) - maximum chromosome length in the fasta file  ]

    main:

    cram_crai = reads_fasta.map { meta, cram, crai, interval, meta_fasta, _fasta, _fai ->
        [
            [
                id: meta.id + "_" + meta_fasta.id,
                sample: meta.id,
                type: meta.datatype,
                fasta_id: meta_fasta.id.tokenize(".")[0..-2].join("."),
            ],
            cram,
            crai,
            interval,
        ]
    }

    // fasta
    fasta = reads_fasta.map { _meta, _cram, _crai, _interval, meta_fasta, fasta, _fai -> [meta_fasta, fasta] }

    // fai
    fai = reads_fasta.map { _meta, _cram, _crai, _interval, meta_fasta, _fasta, fai -> [meta_fasta, fai] }

    // split fasta in compressed format, no gzi index file needed
    gzi = [[], []]
    par_bed = [[], []]

    // call deepvariant
    DEEPVARIANT(cram_crai, fasta, fai, gzi, par_bed)

    // group the vcf files together by sample
    vcf = DEEPVARIANT.out.vcf
        .join(DEEPVARIANT.out.vcf_index)
        .map { meta, vcf, index ->
            [
                [id: [meta.fasta_id, meta.type, meta.sample].join(".")],
                vcf,
                index,
            ]
        }
        .groupTuple()

    // concat vcf files
    BCFTOOLS_CONCAT_VCF(vcf)

    // group the g vcf files together by sample
    g_vcf = DEEPVARIANT.out.gvcf
        .join(DEEPVARIANT.out.gvcf_index)
        .map { meta, gvcf, index ->
            [
                [id: [meta.fasta_id, meta.type, meta.sample].join(".")],
                gvcf,
                index,
            ]
        }
        .groupTuple()

    // concat g vcf files
    BCFTOOLS_CONCAT_GVCF(g_vcf)

    // compress the vcf and gvcf files
    vcf_to_compress = BCFTOOLS_CONCAT_VCF.out.vcf
        .mix(BCFTOOLS_CONCAT_GVCF.out.vcf)
        .combine(max_length)
    ch_compressed_vcf = BGZIPTABIX(vcf_to_compress, [[], [], []]).gz_index

    // generate vcf stats report
    VCF_STATS_REPORT(BCFTOOLS_CONCAT_VCF.out.vcf)

    // generate g vcf stats report
    GVCF_STATS_REPORT(BCFTOOLS_CONCAT_GVCF.out.vcf)

    emit:
    vcf               = BCFTOOLS_CONCAT_VCF.out.vcf // channel: [ val(meta), path(vcf) ]
    gvcf              = BCFTOOLS_CONCAT_GVCF.out.vcf // channel: [ val(meta), path(gvcf) ]
    compressed_vcf    = ch_compressed_vcf // channel: [ val(meta), path(vcf), path(gzi)]
    vcf_csi           = BGZIPTABIX.out.csi // channel: [ val(meta), path(csi)]
    vcf_tbi           = BGZIPTABIX.out.tbi // channel: [ val(meta), path(tbi)]
    vcf_stats_report  = VCF_STATS_REPORT.out.report // channel: [ val(meta), path(report) ]
    gvcf_stats_report = GVCF_STATS_REPORT.out.report // channel: [ val(meta), path(report) ]
}
