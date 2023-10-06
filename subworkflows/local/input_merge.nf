//
// Merge READS(bam or cram files) together by sample name
//

include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge'
include { SAMTOOLS_SORT }  from '../../modules/nf-core/samtools/sort'

workflow INPUT_MERGE {
    take:
    fasta              // file: /path/to/genome.fasta or /path/to/genome.fasta.gz
    fai                // file: /path/to/genome.*.fai
    gzi                // file: /path/to/genome.fasta.gz.gzi or null
    reads              // channel: [ val(meta), data ]
    sort_input         // bollean: true or false

    main:
    ch_versions = Channel.empty()
    
    // sort input reads if asked
    if ( sort_input ) {
      SAMTOOLS_SORT( reads )
      ch_versions = ch_versions.mix ( SAMTOOLS_SORT.out.versions )
      sorted_reads = SAMTOOLS_SORT.out.bam
    } else {     
      sorted_reads = reads
    }

    // group input reads file by sample name
    sorted_reads
     .map{ meta, bam_cram -> [ meta.sample, bam_cram ] }
     .groupTuple()
     .set{ merged_reads } 

    // group input meta data together by sample name as well
    // use the first meta data for the combined reads
    reads
     .map{ meta, bam_cram -> [ meta.sample, meta ] }
     .groupTuple()
     .map { sample, meta_list -> [sample, meta_list[0]] }
     .join( merged_reads )
     .map { sample, meta, bam_cram_list -> [ 
          [ id: ( bam_cram_list.size() == 1 ) ? sample : sample + '_combined',
            type: meta.type 
          ], 
            bam_cram_list 
          ]}
     .set { merged_reads_with_meta }

    // call samtool merge
    ch_fasta = fasta.map { fasta -> [ [ 'id': fasta.baseName ], fasta ] }.first()
    ch_fai = fai.map { fai -> [ [ 'id': fai.baseName ], fai ] }.first()
    ch_gzi = gzi.map { gzi -> [ [ 'id': gzi.baseName ], gzi ] }.first()

    SAMTOOLS_MERGE( merged_reads_with_meta, 
                    ch_fasta,
                    ch_fai,
                    ch_gzi
    )
    ch_versions = ch_versions.mix ( SAMTOOLS_MERGE.out.versions )

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
