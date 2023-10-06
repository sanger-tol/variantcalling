//
// Merge READS(bam or cram files) together by sample name
//

include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge'

workflow INPUT_MERGE {
    take:
    fasta              // file: /path/to/genome.fasta or /path/to/genome.fasta.gz
    fai                // file: /path/to/genome.*.fai
    gzi                // file: /path/to/genome.fasta.gz.gzi or null
    reads              // channel: [ val(meta), data ]

    main:
    // group input reads file by sample name
    reads
     .map{ it -> [ it[0].sample, it[1] ] }
     .groupTuple()
     .set{ merged_reads } 
    
    // group input meta data together by sample name as well
    // use the first meta data for the combined reads
    reads
     .map{ it -> [ it[0].sample, it[0] ] }
     .groupTuple()
     .map { it -> [it[0], it[1][0]] }
     .join( merged_reads )
     .map { it -> [ 
          [ id: ( it[2].size() == 1 ) ? it[1].sample : it[1].sample + '_combined',
            type: it[1].type 
          ], 
            it[2] 
          ]}
     .set { merged_reads_with_meta }

    // call samtool merge
    SAMTOOLS_MERGE( merged_reads_with_meta, 
                    [ [], fasta ],
                    [ [], fai ],
                    [ [], gzi ]
    )

    emit:
    bam      = SAMTOOLS_MERGE.out.bam
    cram     = SAMTOOLS_MERGE.out.cram 
    csi      = SAMTOOLS_MERGE.out.csi
    versions = SAMTOOLS_MERGE.out.versions // channel: [ versions.yml ]

}
