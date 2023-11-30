//
// Merge READS(bam or cram files) together by sample name
//

include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge'
include { SAMTOOLS_SORT }  from '../../modules/nf-core/samtools/sort'

workflow INPUT_MERGE {
    take:
    fasta              // channel: [ val(meta), /path/to/genome.fasta or /path/to/genome.fasta.gz ]
    fai                // channel: [ val(meta), /path/to/genome.*.fai or /path/to/genome.fasta.gz.gzi ]
    reads              // channel: [ val(meta), data ]

    main:
    ch_versions = Channel.empty()

    // group input meta data together by sample name
    reads
     .map{ meta, bam_cram -> [ meta.sample, meta ] }
     .groupTuple()
     .set{ grouped_reads_meta }

    // sort input reads
    SAMTOOLS_SORT( reads )
    ch_versions = ch_versions.mix ( SAMTOOLS_SORT.out.versions )
    sorted_reads = SAMTOOLS_SORT.out.bam

    // group input reads file by sample name
    sorted_reads
     .map{ meta, bam_cram -> [ meta.sample, bam_cram ] }
     .groupTuple()
     .set{ grouped_reads } 

    // join grouped reads and meta
    // use the first meta data for the combined reads
    grouped_reads_meta 
     .map { sample, meta_list -> [sample, meta_list[0]] }
     .join( grouped_reads )
     .map { sample, meta, bam_cram_list -> [ 
          [ id: ( bam_cram_list.size() == 1 ) ? sample : sample + '_combined',
            datatype: meta.datatype 
          ], 
            bam_cram_list 
          ]}
     .set { grouped_reads_with_meta }

    // call samtool merge
    SAMTOOLS_MERGE( grouped_reads_with_meta, 
                    fasta,
                    fai
    )
    ch_versions = ch_versions.mix ( SAMTOOLS_MERGE.out.versions )

    // concat merged bam or cram together along with their index file
    SAMTOOLS_MERGE.out.bam
      .join(SAMTOOLS_MERGE.out.csi)
      .concat(
        SAMTOOLS_MERGE.out.cram
         .join(SAMTOOLS_MERGE.out.crai)
      )
    .set{ indexed_merged_reads };

    emit:
    indexed_merged_reads = indexed_merged_reads
    versions = ch_versions // channel: [ versions.yml ]

}
