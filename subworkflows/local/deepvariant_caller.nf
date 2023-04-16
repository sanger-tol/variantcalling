//
// call deepvariant
//

include { DEEPVARIANT }         from '../../modules/nf-core/deepvariant/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_VCF   }  from '../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_CONCAT as BCFTOOLS_CONCAT_GVCF }   from '../../modules/nf-core/bcftools/concat/main'

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

    // group the vcf files together by sample
    DEEPVARIANT.out.vcf
     .map{ [ it[0].sample, it[1] ] }
     .groupTuple()
     .map { [ [id: it[0]], it[1], [] ]}
     .set{ vcf }
    
    // catcat vcf files
    BCFTOOLS_CONCAT_VCF( vcf )
    ch_versions = ch_versions.mix ( BCFTOOLS_CONCAT_VCF.out.versions.first() )

    // group the g vcf files together by sample
    DEEPVARIANT.out.gvcf
     .map{ [ it[0].sample, it[1] ] }
     .groupTuple()
     .map { [ [ id: it[0] + '.g' ], it[1], [] ]}
     .set{ g_vcf }
    
    // catcat g vcf files
    BCFTOOLS_CONCAT_GVCF( g_vcf )
    ch_versions = ch_versions.mix ( BCFTOOLS_CONCAT_GVCF.out.versions.first() )

    emit:
    vcf      = BCFTOOLS_CONCAT_VCF.out.vcf         // /path/to/vcf
    gvcf     = BCFTOOLS_CONCAT_GVCF.out.vcf        // /path/to/gvcf
    versions = ch_versions                         // channel: [ versions.yml ]
}
