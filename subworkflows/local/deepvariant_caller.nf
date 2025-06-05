//
// Call variants with Deepvariant
//

include { DEEPVARIANT_RUNDEEPVARIANT as DEEPVARIANT }   from '../../modules/nf-core/deepvariant/rundeepvariant/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF    }   from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_GVCF   }   from '../../modules/nf-core/bcftools/concat/main'
include { DEEPVARIANT_VCFSTATSREPORT as VCF_STATS_REPORT }   from '../../modules/nf-core/deepvariant/vcfstatsreport/main'
include { DEEPVARIANT_VCFSTATSREPORT as GVCF_STATS_REPORT }   from '../../modules/nf-core/deepvariant/vcfstatsreport/main'

workflow DEEPVARIANT_CALLER {
    take:
    reads_fasta    // [ val(meta), cram, crai, interval, fasta_file_name, fasta, fai ]

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

    // generate vcf stats report
    VCF_STATS_REPORT ( BCFTOOLS_CONCAT_VCF.out.vcf )
    ch_versions = ch_versions.mix ( VCF_STATS_REPORT.out.versions.first() )

    // generate g vcf stats report
    GVCF_STATS_REPORT ( BCFTOOLS_CONCAT_GVCF.out.vcf )
    ch_versions = ch_versions.mix ( GVCF_STATS_REPORT.out.versions.first() )

    emit:
    vcf      = BCFTOOLS_CONCAT_VCF.out.vcf         // channel: [ val(meta), path(vcf) ]
    gvcf     = BCFTOOLS_CONCAT_GVCF.out.vcf        // channel: [ val(meta), path(gvcf) ]
    versions = ch_versions                         // channel: [ versions.yml ]
}
