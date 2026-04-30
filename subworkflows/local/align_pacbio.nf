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
    fastq_filter = FILTER_PACBIO.out.fastq
        .combine(fasta)
        .map { meta, fastq, meta_fasta, fasta, _fai -> [meta + ['fasta_id': meta_fasta.id, 'id':"${meta_fasta.id}.${meta.datatype}.${meta.id}"], fastq]}
    MINIMAP2_ALIGN(fastq_filter, fasta.map { meta, fa, _fai -> [meta, fa] }, true, false, false, false)

    // Collect all alignment output by sample name
    ch_bams = MINIMAP2_ALIGN.out.bam
        .map { meta, bam -> [['id': meta.specimen, 'datatype': meta.datatype], [['id':meta.id, 'specimen': meta.specimen, 'datatype': meta.datatype, 'sample': meta.sample, 'run': meta.run, 'fasta_id': meta.fasta_id], bam]] }
        .groupTuple(by: [0])
        .branch { meta, bams ->
            to_merge: bams.size() > 1
            no_merge: true
        }

    ch_bams_no_merge = ch_bams.no_merge
        .map { meta, bams -> [ bams[0][0], bams[0][1], [] ] }

    ch_bams_to_merge = ch_bams.to_merge
        .map { meta, orig_id_bams ->
            def meta_bam = orig_id_bams[0][0]
            def meta_bam_new = meta_bam + ['sample': "${meta_bam.specimen}/merge", 'id': "${meta_bam.fasta_id}.${meta_bam.datatype}.${meta_bam.specimen}.merge", 'run': "merge"]
            def bams = orig_id_bams
                .sort { a, b -> a[0].id <=> b[0].id} // sort by id to ensure consistent order
                .collect { id_bam -> id_bam[1] }
            [meta_bam_new, bams, []]
        }


    // Merge
    SAMTOOLS_MERGE(ch_bams_to_merge, [[], [], [], []])


    // Convert merged BAM to CRAM and calculate indices and statistics
    ch_sort = SAMTOOLS_MERGE.out.bam.map { meta, bam -> [meta, bam, []] }.mix(ch_bams_no_merge)
    CONVERT_STATS(ch_sort, fasta)

    emit:
    cram     = CONVERT_STATS.out.cram // channel: [ val(meta), /path/to/cram ]
    crai     = CONVERT_STATS.out.crai // channel: [ val(meta), /path/to/crai ]
    stats    = CONVERT_STATS.out.stats // channel: [ val(meta), /path/to/stats ]
    idxstats = CONVERT_STATS.out.idxstats // channel: [ val(meta), /path/to/idxstats ]
    flagstat = CONVERT_STATS.out.flagstat // channel: [ val(meta), /path/to/flagstat ]
}
