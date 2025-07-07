//
// Run Himut
//

include { HIMUT                    }   from '../../modules/local/himut.nf'
include { TABIX_BGZIP as BGZIP     }   from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as TABIX_CSI }   from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TBI }   from '../../modules/nf-core/tabix/tabix/main'

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
    HIMUT (
        ch_genome,
        ch_region_list,
        ch_aligned_reads,
        ch_aligned_reads_index,
        vcf_input,
        vcf_index
    )
    ch_versions = ch_versions.mix ( HIMUT.out.versions.first() )

    // compress the himut vcf files
    ch_compressed_himut_vcf = BGZIP ( HIMUT.out.vcf_output ).output
    ch_versions             = ch_versions.mix ( BGZIP.out.versions.first() )



    // index the compressed himut vcf files
    ch_compressed_himut_vcf
        .combine(max_length)
        .map { meta_vcf, vcf, meta -> [ meta_vcf + meta, vcf ] }
        .branch { meta, vcf ->
        tbi_and_csi: meta.max_length < 2**29
        only_csi:    meta.max_length < 2**32
        }
        .set { himut_tabix_selector }

    // do the indexing on the compatible gvcf files
    ch_himut_vcf_csi = TABIX_CSI ( himut_tabix_selector.tbi_and_csi.mix(himut_tabix_selector.only_csi) ).csi
    ch_versions        = ch_versions.mix ( TABIX_CSI.out.versions.first() )
    ch_himut_vcf_tbi = TABIX_TBI ( himut_tabix_selector.tbi_and_csi ).tbi
    ch_versions        = ch_versions.mix ( TABIX_TBI.out.versions.first() )



    emit:
    himut_vcf = HIMUT.out.vcf_output
    compressed_himut_vcf = ch_compressed_himut_vcf // channel: [ val(meta), path(output)]
    himut_vcf_csi  = ch_himut_vcf_csi              // channel: [ val(meta), path(csi)]
    himut_vcf_tbi  = ch_himut_vcf_tbi              // channel: [ val(meta), path(tbi)]
    versions = ch_versions
}



