//
// Check input samplesheet and get read channels
//

include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'

workflow INPUT_FILTER_SPLIT {

    take:
    fasta    // path to reference fasta file, either compressed or uncompressed
    fai      // path to compressed or uncompressed reference fasta index file
    gzi      // path to compressed fasta index file

    main:
    ch_versions = Channel.empty()

    // split the fasta file into files with one sequence each
    Channel
     .fromPath(fasta)
     .splitFasta(file:true)
     .map{ [ [id: 'splitting_fasta'],  it ] }
     .set{ splitted_fasta }

    // index splitted fasta files
    SAMTOOLS_FAIDX ( splitted_fasta  )
    ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )
   
    splitted_fasta.join( SAMTOOLS_FAIDX.out.fai ).view()

    
    emit:        
    versions   = ch_versions                       // channel: [ versions.yml ]

}
