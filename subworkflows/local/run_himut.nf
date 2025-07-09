//
// Run Himut
//

include { HIMUT                    }   from '../../modules/local/himut.nf'
include { TABIX_BGZIP as BGZIP     }   from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as TABIX_CSI }   from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TBI }   from '../../modules/nf-core/tabix/tabix/main'

workflow RUN_HIMUT {
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

    // Transfer .cram inout to .bam
    ch_aligned_reads = ch_aligned_reads.map { meta, bam ->
        if (bam.endsWith('.cram')) {
            def bam_index = ch_aligned_reads_index.find { it.meta.id == meta.id }
            def bam_out = "${bam}.bam"
            def cmd = "samtools view -b -T ${bam_index} ${bam} > ${bam_out}"
            log.info "Converting ${bam} to ${bam_out}"
            sh(cmd)
            return [meta, bam_out]
        } else {
            return [meta, bam]
        }
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



