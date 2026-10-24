//
// Filter PacBio reads
// Original protocol is a modified version by Shane of the original program, HiFiAdapterFilt
//

include { SAMTOOLS_VIEW as SAMTOOLS_CONVERT } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_COLLATE                  } from '../../modules/nf-core/samtools/collate/main'
include { SAMTOOLS_FASTA                    } from '../../modules/nf-core/samtools/fasta/main'
include { GUNZIP                            } from '../../modules/nf-core/gunzip/main'
include { BLAST_BLASTN                      } from '../../modules/nf-core/blast/blastn/main'
include { PACBIO_FILTER                     } from '../../modules/local/pacbio_filter'
include { SAMTOOLS_VIEW as SAMTOOLS_FILTER  } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ                    } from '../../modules/nf-core/samtools/fastq/main'


workflow FILTER_PACBIO {
    take:
    reads // channel: [ val(meta), /path/to/datafile ]
    db // channel: /path/to/vector_db

    main:
    // Convert from PacBio BAM to Samtools BAM
    ch_pacbio = reads.map { meta, bam -> [meta, bam, []] }

    SAMTOOLS_CONVERT(ch_pacbio, [[], [], []], [[], []], [[], []], [])


    // Collate BAM file to create interleaved FASTA
    SAMTOOLS_COLLATE(SAMTOOLS_CONVERT.out.bam, [[], [], []])


    // Convert BAM to FASTA
    SAMTOOLS_FASTA(SAMTOOLS_COLLATE.out.bam, true)


    // Gunzip FASTA file to BLAST
    GUNZIP(SAMTOOLS_FASTA.out.other)


    // Nucleotide BLAST
    ch_db = db.map { path -> [[], path] }
    BLAST_BLASTN(GUNZIP.out.gunzip, ch_db, [], [], [])


    // Filter BLAST output
    PACBIO_FILTER(BLAST_BLASTN.out.txt)


    // Create filtered BAM file
    ch_reads_and_list = SAMTOOLS_CONVERT.out.bam
        .join(SAMTOOLS_CONVERT.out.csi)
        .join(PACBIO_FILTER.out.list)

    ch_reads = ch_reads_and_list.map { meta, bam, csi, _list -> [meta, bam, csi] }

    ch_lists = ch_reads_and_list.map { meta, _bam, _csi, list -> [meta, list] }

    SAMTOOLS_FILTER(ch_reads, [[], [], []], ch_lists, [[], []], [])


    // Convert BAM to FASTQ
    SAMTOOLS_FASTQ(SAMTOOLS_FILTER.out.unselected, true)

    emit:
    fastq    = SAMTOOLS_FASTQ.out.other // channel: [ meta, /path/to/fastq ]
}
