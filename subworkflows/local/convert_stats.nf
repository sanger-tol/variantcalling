//
// Convert BAM to CRAM, create index and calculate statistics
//

include { SAMTOOLS_VIEW     } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_STATS    } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/samtools/idxstats/main'


workflow CONVERT_STATS {
    take:
    bam // channel: [ val(meta), /path/to/bam, /path/to/bai]
    fasta // channel: [ val(meta), /path/to/fasta[.gz], /path/to/fai ]

    main:
    // Convert BAM to CRAM
    SAMTOOLS_VIEW(bam, fasta, [[],[]], [[],[]], [])


    // Combine CRAM and CRAI into one channel
    ch_cram_crai = SAMTOOLS_VIEW.out.cram.join(SAMTOOLS_VIEW.out.crai)


    // Calculate statistics
    SAMTOOLS_STATS(ch_cram_crai, fasta)


    // Calculate statistics based on flag values
    SAMTOOLS_FLAGSTAT(ch_cram_crai)


    // Calculate index statistics
    SAMTOOLS_IDXSTATS(ch_cram_crai)

    emit:
    cram     = SAMTOOLS_VIEW.out.cram // channel: [ val(meta), /path/to/cram ]
    crai     = SAMTOOLS_VIEW.out.crai // channel: [ val(meta), /path/to/crai ]
    stats    = SAMTOOLS_STATS.out.stats // channel: [ val(meta), /path/to/stats ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), /path/to/idxstats ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), /path/to/flagstat ]
}
