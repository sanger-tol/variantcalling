//
// Align PacBio read files against the genome
//

include { FILTER_PACBIO  } from '../../subworkflows/local/filter_pacbio'
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main'
include { CONVERT_STATS  } from '../../subworkflows/local/convert_stats'


workflow ALIGN_PACBIO {
    take:
    fasta // channel: [ val(meta), /path/to/fasta[.gz], /path/to/fai ]
    reads // channel: [ val(meta), /path/to/datafile ]
    db // channel: /path/to/vector_db

    main:
    // Filter BAM and output as FASTQ
    FILTER_PACBIO(reads, db)


    // Align Fastq to Genome
    MINIMAP2_ALIGN(FILTER_PACBIO.out.fastq, fasta.map { meta, fa, _fai -> [meta, fa] }, true, false, false, false)


    // Collect all alignment output by sample name
    ch_bams = MINIMAP2_ALIGN.out.bam
        .map { meta, bam -> [['id': meta.sample, 'datatype': meta.datatype, 'sample': meta.sample], [meta.id, bam]] }
        .groupTuple(by: [0])
        .map { meta, orig_id_bams ->
            def bams = orig_id_bams
                .sort { a, b -> a[0] <=> b[0]} // sort by id to ensure consistent order
                .collect { id_bam -> id_bam[1] }
            [meta, bams, []]
        }


    // Merge
    SAMTOOLS_MERGE(ch_bams, [[], [], [], []])


    // Convert merged BAM to CRAM and calculate indices and statistics
    ch_sort = SAMTOOLS_MERGE.out.bam.map { meta, bam -> [meta, bam, []] }
    CONVERT_STATS(ch_sort, fasta)

    emit:
    cram     = CONVERT_STATS.out.cram // channel: [ val(meta), /path/to/cram ]
    crai     = CONVERT_STATS.out.crai // channel: [ val(meta), /path/to/crai ]
    stats    = CONVERT_STATS.out.stats // channel: [ val(meta), /path/to/stats ]
    idxstats = CONVERT_STATS.out.idxstats // channel: [ val(meta), /path/to/idxstats ]
    flagstat = CONVERT_STATS.out.flagstat // channel: [ val(meta), /path/to/flagstat ]
}
