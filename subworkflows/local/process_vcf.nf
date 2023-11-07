//
// Call vcftools to process VCF files
//

include { VCFTOOLS as VCFTOOLS_SITE_PI   }   from '../../modules/nf-core/vcftools/main'
include { VCFTOOLS as VCFTOOLS_HET       }   from '../../modules/nf-core/vcftools/main'

workflow PROCESS_VCF {
    take:
    vcf    // [ val(meta), vcf ]

    main:
    ch_versions = Channel.empty()

    // call vcftools for per site nucleotide diversity
    VCFTOOLS_SITE_PI(
      vcf, [], []
    )
    ch_versions = ch_versions.mix( VCFTOOLS_SITE_PI.out.versions )

    // call vcftools to calculate for heterozygosity
    VCFTOOLS_HET(
      vcf, [], []
    )
    ch_versions = ch_versions.mix( VCFTOOLS_HET.out.versions )

    emit:
    versions       = ch_versions                      // channel: [ versions.yml ]
    stite_pi       = VCFTOOLS_SITE_PI.out.sites_pi    // [ meta, site_pi ]
    heterozygosity = VCFTOOLS_HET.out.heterozygosity  // [ meta,  heterozygosity ]

}
