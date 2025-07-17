//
// Convert to BAM, create index and calculate statistics
//

include { SAMTOOLS_VIEW     } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_STATS    } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/samtools/idxstats/main'


workflow CONVERT_STATS {
    take:
    bam       // channel: [ val(meta), /path/to/bam, /path/to/bai]
    fasta     // channel: [ val(meta), /path/to/fasta ]


    main:
    ch_versions = Channel.empty()


    // Convert input to BAM
    SAMTOOLS_VIEW ( bam, fasta, [ ] )
    ch_versions = ch_versions.mix ( SAMTOOLS_VIEW.out.versions.first() )


    // Combine BAM and BAI into one channel
    SAMTOOLS_VIEW.out.bam
        .join ( SAMTOOLS_VIEW.out.bai.mix(SAMTOOLS_VIEW.out.csi) )
        .set { ch_bam_bai }


    // Calculate statistics
    SAMTOOLS_STATS ( ch_bam_bai, fasta )
    ch_versions = ch_versions.mix ( SAMTOOLS_STATS.out.versions.first() )


    // Calculate statistics based on flag values
    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions = ch_versions.mix ( SAMTOOLS_FLAGSTAT.out.versions.first() )


    // Calculate index statistics
    SAMTOOLS_IDXSTATS ( ch_bam_bai )
    ch_versions = ch_versions.mix ( SAMTOOLS_IDXSTATS.out.versions.first() )


    emit:
    bam      = SAMTOOLS_VIEW.out.bam             // channel: [ val(meta), /path/to/bam      ]
    csi      = SAMTOOLS_VIEW.out.csi             // channel: [ val(meta), /path/to/csi      ][ SAMTOOLS_VIEW.out.bai in empty ]
    stats    = SAMTOOLS_STATS.out.stats          // channel: [ val(meta), /path/to/stats    ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat    // channel: [ val(meta), /path/to/idxstats ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats    // channel: [ val(meta), /path/to/flagstat ]
    versions = ch_versions                       // channel: [ versions.yml ]
}
