//
// call deepvariant
//

include { DEEPVARIANT }         from '../../modules/nf-core/deepvariant/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF   }   from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_GVCF  }   from '../../modules/nf-core/bcftools/concat/main'

workflow DEEPVARIANT_CALLER {

    take:
    reads_fasta    // [ val(meta), cram, crai, interval, fasta_file_name, fasta, fai ]

    main:
    ch_versions = Channel.empty()

    reads_fasta.map{ meta, cram, crai, interval, fasta_file_name, fasta, fai ->
                    [ [ id: meta.id + "_" + fasta_file_name, sample: meta.id, type: meta.type ], 
                       cram,
                       crai, 
                       interval 
                    ] }
               .set{ cram_crai } // [ val(meta), cram, crai, interval ]

    // fasta
    reads_fasta.map{ meta, cram, crai, interval, fasta_file_name, fasta, fai -> [ fasta ] }
            .set{ fasta }

    // fai
    reads_fasta.map{ meta, cram, crai, interval, fasta_file_name, fasta, fai -> [ fai ] }
            .set{ fai }

    // split fasta in compressed format, no gzi index file needed
    gzi = []

    // call deepvariant
    DEEPVARIANT( cram_crai, fasta, fai, gzi )
    ch_versions = ch_versions.mix ( DEEPVARIANT.out.versions.first() )

    // group the vcf files together by sample
    DEEPVARIANT.out.vcf
     .map{ meta, vcf -> [ meta.sample, vcf ] }
     .groupTuple()
     .map { sample, vcf -> [ [id: sample], vcf, [] ]}
     .set{ vcf }
    
    // catcat vcf files
    BCFTOOLS_CONCAT_VCF( vcf )
    ch_versions = ch_versions.mix ( BCFTOOLS_CONCAT_VCF.out.versions.first() )

    // group the g vcf files together by sample
    DEEPVARIANT.out.gvcf
     .map{ meta, gvcf -> [ meta.sample, gvcf ] }
     .groupTuple()
     .map { sample, gvcf ->  [ [ id: sample + '.g' ], gvcf, [] ]}
     .set{ g_vcf }
    
    // catcat g vcf files
    BCFTOOLS_CONCAT_GVCF( g_vcf )
    ch_versions = ch_versions.mix ( BCFTOOLS_CONCAT_GVCF.out.versions.first() )

    emit:
    vcf      = BCFTOOLS_CONCAT_VCF.out.vcf         // /path/to/vcf
    gvcf     = BCFTOOLS_CONCAT_GVCF.out.vcf        // /path/to/gvcf
    versions = ch_versions                         // channel: [ versions.yml ]
}
