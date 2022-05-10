//
// Sub-workflow to perform variant calling on illumina data
//

include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main'

workflow ILLUMINA_VC {
    take:
    reference_fasta     // channel: path/to/reference_fasta
    mapped_cram         // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()

    //
    // gatk4 Create Dictionary for the Reference FASTA
    //
    GATK4_CREATESEQUENCEDICTIONARY ( reference_fasta )
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions.first())

}