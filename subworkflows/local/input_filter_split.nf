//
// Split input fasta file by sequence and filter the input reads
//

include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_VIEW  } from '../../modules/nf-core/samtools/view/main'
include { CAT_CAT        } from '../../modules/nf-core/cat/cat/main'

workflow INPUT_FILTER_SPLIT {
    take:
    fasta              // file: /path/to/genome.fasta or /path/to/genome.fasta.gz
    fai                // file: /path/to/genome.*.fai
    gzi                // file: /path/to/genome.fasta.gz.gzi or null
    reads              // [ val(meta), data, index ]
    interval           // file: /path/to/intervals.bed
    split_fasta_cutoff // val(min_file_size)

    main:
    ch_versions = Channel.empty()

    // split the fasta file into files with one sequence each, group them by file size
    Channel
     .fromPath ( fasta )
     .splitFasta ( file:true )
     .branch {
        small: it.size() < split_fasta_cutoff
        large: it.size() >= split_fasta_cutoff
     }
     .set { branched_fasta_files }
     
    // check the large split fasta files
    branched_fasta_files.large
     .map { large_file -> [ [ id: large_file.baseName ], large_file ] }
     .set { ch_large_files }
    
    // check all the small split fasta files
    branched_fasta_files.small
     .collect()
     .map { small_files -> [ [ id : small_files[0].baseName.substring(0, small_files[0].baseName.lastIndexOf('.') ) + '.small' ], small_files ] }
     .set { ch_samll_files }
    
    // merge all small split fasta files together
    CAT_CAT ( ch_samll_files )
    ch_versions = ch_versions.mix ( CAT_CAT.out.versions )
    
    // concat large and merged samll fasta files together
    Channel.empty()
     .concat ( CAT_CAT.out.file_out, ch_large_files )
     .set { split_fasta } 

    // index split fasta files
    SAMTOOLS_FAIDX ( split_fasta,  [[], []])
    ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions.first() )
    
    // join fasta with corresponding fai file
    split_fasta
     .map { meta, fasta -> [ fasta.baseName, fasta ] }
     .join ( 
        SAMTOOLS_FAIDX.out.fai
         .map { mata, fai -> [ fai.baseName - '.fasta', fai ] } 
      )
     .set { fasta_fai }

    // filter reads
    SAMTOOLS_VIEW ( reads, [ [], fasta ], [] )
    ch_versions = ch_versions.mix ( SAMTOOLS_VIEW.out.versions.first() )
    
    // combine reads with splitted references
    SAMTOOLS_VIEW.out.cram
     .join ( SAMTOOLS_VIEW.out.crai )
     .map { filtered_reads -> filtered_reads + [interval ?: []] }
     .combine ( fasta_fai )
     .set { cram_crai_fasta_fai }

    emit:
    reads_fasta    = cram_crai_fasta_fai  // channel: [ val(meta), cram, crai, interval, fasta_file_name, fasta, fai ]
    versions       = ch_versions          // channel: [ versions.yml ]
}
