//
// Split input fasta file by sequence and filter the input reads
//

include { GUNZIP         } from '../../modules/nf-core/gunzip/main'
include { SEQKIT_SPLIT2  } from '../../modules/nf-core/seqkit/split2/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view/main'

workflow INPUT_FILTER_SPLIT {
    take:
    fasta // [ val(meta, /path/to/fasta[.gz], /path/to/fai) ]
    reads // [ val(meta), data, index ]
    interval // file: /path/to/intervals.bed

    main:
    //
    // MODULE: Unzip the fasta if zipped
    //
    ch_fasta = fasta.branch { meta, fa, _fai ->
        gzipped: fa.name.endsWith('.gz')
            [meta, fa]
        unzipped: true
            [meta, fa]
    }

    GUNZIP(
        ch_fasta.gzipped
    )

    ch_fasta_to_split = GUNZIP.out.gunzip.mix(ch_fasta.unzipped)

    //
    // MODULE: Split the Fasta file in chunks
    //
    SEQKIT_SPLIT2(ch_fasta_to_split)

    // Add pertinent meta maps to the chunks
    ch_split_fastas = SEQKIT_SPLIT2.out.reads
        .map { _meta, fastas -> fastas }
        .flatten()
        .map { fa -> [[id: fa.baseName, total_length: fa.size()], fa, []] }

    //
    // MODULE: Index the chunks
    //
    SAMTOOLS_FAIDX(ch_split_fastas, [[], []])

    // join fasta with corresponding fai file
    fasta_fai = ch_split_fastas.join(SAMTOOLS_FAIDX.out.fai)

    // filter reads
    ch_fasta = fasta.map { _meta, fa, fai -> [['id': fa.baseName], fa, fai] }.first()

    SAMTOOLS_VIEW(reads, ch_fasta, [[], []], [[], []], [])

    // combine reads with splitted references
    cram_crai_fasta_fai = SAMTOOLS_VIEW.out.cram
        .join(SAMTOOLS_VIEW.out.crai)
        .combine(interval.ifEmpty([[]]))
        .combine(fasta_fai)

    emit:
    reads_fasta = cram_crai_fasta_fai // channel: [ val(meta), cram, crai, interval, val(meta_fasta), fasta, fai ]
}
