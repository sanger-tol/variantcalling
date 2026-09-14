//
// Split input fasta file by sequence and filter the input reads
//

include { SEQKIT_SPLIT2  } from '../../modules/nf-core/seqkit/split2/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view/main'

workflow INPUT_FILTER_SPLIT {
    take:
    fasta // [ val(meta, /path/to/fasta[.gz], /path/to/fai) ]
    reads // [ val(meta), data, index ]
    intervals // file: /path/to/intervals.bed

    main:
    //
    // MODULE: Split the Fasta file in chunks
    //
    ch_fasta_for_split = fasta.map { meta, fa, fai -> [meta, fa] }
    SEQKIT_SPLIT2(ch_fasta_for_split)

    // Add pertinent meta maps to the chunks
    ch_split_fastas = SEQKIT_SPLIT2.out.reads
        .map { _meta, fastas -> fastas }
        .flatten()
        .map { fa -> [[id: fa.baseName, total_length: fa.size()], fa, []] }

    //
    // MODULE: Index the chunks
    //
    SAMTOOLS_FAIDX(ch_split_fastas, false)

    // join fasta with corresponding fai file
    fasta_fai = ch_split_fastas.join(SAMTOOLS_FAIDX.out.fai)

    // filter reads
    ch_fasta = fasta.map { _meta, fa, fai -> [['id': fa.baseName], fa, fai] }.first()

    SAMTOOLS_VIEW(reads, ch_fasta, [[], []], [[], []], [])

    // combine reads with splitted references
    cram_crai_fasta_fai = SAMTOOLS_VIEW.out.cram
        .join(SAMTOOLS_VIEW.out.crai)
        .combine(intervals.ifEmpty([[]]))
        .combine(fasta_fai)

    emit:
    reads_fasta = cram_crai_fasta_fai // channel: [ val(meta), cram, crai, intervals, val(meta_fasta), fasta, [], fai ]
}
