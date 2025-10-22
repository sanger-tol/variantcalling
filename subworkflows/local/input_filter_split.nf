//
// Split input fasta file by sequence and filter the input reads
//

include { GUNZIP         } from '../../modules/nf-core/gunzip/main'
include { SEQKIT_SPLIT2  } from '../../modules/nf-core/seqkit/split2/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view/main'

workflow INPUT_FILTER_SPLIT {
    take:
    fasta              // [ val(meta, /path/to/genome.fasta[.gz] ]
    reads              // [ val(meta), data, index ]
    interval           // file: /path/to/intervals.bed

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Unzip the fasta if zipped
    //
    fasta
    | branch { meta, fa ->
        gzipped: fa.name.endsWith('.gz')
        unzipped: true
    }
    | set { ch_fasta }

    GUNZIP (
        ch_fasta.gzipped
    )
    ch_versions  = ch_versions.mix ( GUNZIP.out.versions )

    GUNZIP.out.gunzip
    | mix ( ch_fasta.unzipped )
    | set { ch_fasta_to_split }

    //
    // MODULE: Split the Fasta file in chunks
    //
    SEQKIT_SPLIT2 ( ch_fasta_to_split )
    ch_versions = ch_versions.mix ( SEQKIT_SPLIT2.out.versions )

    // Add pertinent meta maps to the chunks
    SEQKIT_SPLIT2.out.reads
    | map { meta, fastas -> fastas }
    | flatten
    | map { fa -> [ [id: fa.baseName, total_length: fa.size()], fa ] }
    | set { ch_split_fastas }

    //
    // MODULE: Index the chunks
    //
    SAMTOOLS_FAIDX ( ch_split_fastas,  [[], []])
    ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions.first() )

    // join fasta with corresponding fai file
    ch_split_fastas
    | join ( SAMTOOLS_FAIDX.out.fai )
    | set { fasta_fai }

    //
    // MODULE: filter the reads
    //
    SAMTOOLS_VIEW ( reads, fasta, [] )
    ch_versions = ch_versions.mix ( SAMTOOLS_VIEW.out.versions.first() )

    // combine reads with splitted references
    SAMTOOLS_VIEW.out.cram
    | join ( SAMTOOLS_VIEW.out.crai )
    | combine(interval.ifEmpty([[]]))
    | combine ( fasta_fai )
    | set { cram_crai_fasta_fai }

    emit:
    reads_fasta    = cram_crai_fasta_fai  // channel: [ val(meta), cram, crai, interval, val(meta_fasta), fasta, fai ]
    versions       = ch_versions          // channel: [ versions.yml ]
}
