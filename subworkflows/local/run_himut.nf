//
// Run Himut
//

include { HIMUT }   from '../../modules/local/himut.nf'

workflow HIMUT {
    take:
    ch_genome,              // genome fasta
    ch_genome_index_fai,    // genome index
    ch_region_list,         // region list, generated from .fai file
    ch_aligned_reads,       // bam
    ch_aligned_reads_index, // bam index
    vcf_input,              // vcf input
    vcf_index               // vcf index

    main:
    ch_versions = Channel.empty()

    // If region list is not provided, generate region_list from .fai file
    // Do not include sex chromosomes
    if (params.region_list) {
        ch_region_list = Channel.fromPath(params.region_list)
    } else {
        def make_region_list = { meta, fai ->
            def chr_list = []
            def chr_length = 0
            def chr_name = ''
            fai.splitEachLine('\t') { line ->
                // How to select out sex chromosomes?
                // if (line[0].startsWith('chr') && !line[0].contains('M') && !line[0].contains('Y')) {
                    chr_name = line[0]
                    chr_length = line[1] as int
                    if (chr_length > 0) {
                        chr_list << chr_name
                    }
            }
        }
        return [meta, chr_list]

        ch_genome_index_fai
            .map { meta, fai -> make_region_list(meta, fai) }
            .set { ch_region_list }
    }

    // Run HIMUT
    HIMUT ( ch_genome, ch_region_list, ch_aligned_reads, ch_aligned_reads_index, vcf_input, vcf_index )
    ch_versions = ch_versions.mix ( HIMUT.out.versions.first() )

    emit:
    himut_vcf = HIMUT.out.vcf_output
    versions = ch_versions
}



