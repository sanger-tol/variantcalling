//
// Merge READS(bam or cram files) together by sample name
//

include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge'
include { SAMTOOLS_SORT  } from '../../modules/nf-core/samtools/sort'

workflow INPUT_MERGE {
    take:
    fasta // channel: [ val(meta), /path/to/fasta[.gz], /path/to/[fai,gzi]]
    reads // channel: [ val(meta), data ]

    main:
    // Add fasta id to the reads meta
    reads = reads.combine(fasta).map { meta, reads, meta_fasta, _fasta, _fai -> [meta + ['fasta_id': meta_fasta.id, 'id':"${meta_fasta.id}.${meta.datatype}.${meta.id}"], reads] }

    // sort input reads
    SAMTOOLS_SORT(reads, fasta, [])
    sorted_reads = SAMTOOLS_SORT.out.bam

    grouped_reads_meta = sorted_reads
        .map { meta, reads -> [meta.specimen, [meta, reads]] }
        .groupTuple()
        .branch { _specimen, meta_reads ->
            to_merge: meta_reads.size() > 1
            no_merge: true
        }

    ch_reads_no_merge = grouped_reads_meta.no_merge.map { _meta, reads -> [ reads[0][0], reads[0][1], [] ] }
    ch_reads_to_merge = grouped_reads_meta.to_merge
        .map { _meta, orig_id_reads ->
            def meta_read = orig_id_reads[0][0]
            def runs = orig_id_reads.collect { id_read -> id_read[0].run }
            def meta_read_new = meta_read + ['sample': "${meta_read.specimen}/${params.merge_output}",
                                            'id': "${meta_read.fasta_id}.${meta_read.datatype}.${meta_read.specimen}.${params.merge_output}",
                                            'run': "merge",
                                            'merge_source': runs.sort().join("\n"),
                                            'basename': meta_read.basename.replaceAll(meta_read.run, params.merge_output) ]
            def reads = orig_id_reads
                .sort { a, b -> a[0].id <=> b[0].id} // sort by id to ensure consistent order
                .collect { id_read -> id_read[1] }
            [meta_read_new, reads, []]
        }

    // call samtool merge
    SAMTOOLS_MERGE(
        ch_reads_to_merge,
        fasta.map { meta, fa, fai -> [meta, fa, fai, []] },
    )

    // concat merged bam or cram together along with their index file
    merged_reads = SAMTOOLS_MERGE.out.bam
        .concat(SAMTOOLS_MERGE.out.cram)
        .join(SAMTOOLS_MERGE.out.index)
        .mix(ch_reads_no_merge)

    emit:
    merged_reads = merged_reads
}
