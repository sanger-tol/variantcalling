//
// Call variants with Deepvariant
//

include { DEEPVARIANT_RUNDEEPVARIANT as DEEPVARIANT }   from '../../modules/nf-core/deepvariant/deepvariant_rundeepvariant/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF    }   from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_GVCF   }   from '../../modules/nf-core/bcftools/concat/main'

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
                             [ [ id: meta.id + "_" + fasta_file_name, sample: meta.id, type: meta.datatype ],
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
     .map { meta, vcf -> [
            [ id: meta.fasta_file_name.tokenize(".")[0..-2].join(".")
                  + "." + meta.type
                  + "." + meta.sample
            ],
            vcf
          ] }
     .groupTuple()
     .map { meta, vcf -> [ meta, vcf, [] ] }
     .set { vcf }

    // catcat vcf files
    BCFTOOLS_CONCAT_VCF ( vcf )
    ch_versions = ch_versions.mix ( BCFTOOLS_CONCAT_VCF.out.versions.first() )

    // group the g vcf files together by sample
    DEEPVARIANT.out.gvcf
     .map { meta, gvcf -> [
            [ id: meta.fasta_file_name.tokenize(".")[0..-2].join(".")
                  + "." + meta.type
                  + "." + meta.sample
            ],
            gvcf
          ] }
     .groupTuple()
     .map { meta, gvcf -> [ meta, gvcf, [] ] }
     .set { g_vcf }

    // catcat g vcf files
    BCFTOOLS_CONCAT_GVCF ( g_vcf )
    ch_versions = ch_versions.mix ( BCFTOOLS_CONCAT_GVCF.out.versions.first() )

    emit:
    vcf      = BCFTOOLS_CONCAT_VCF.out.vcf         // channel: [ val(meta), path(vcf) ]
    gvcf     = BCFTOOLS_CONCAT_GVCF.out.vcf        // channel: [ val(meta), path(gvcf) ]
    versions = ch_versions                         // channel: [ versions.yml ]
}
