

params.fasta = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.fasta"
params.fasta_index = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.fasta.fai"
params.assembly_report = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1_assembly_report.txt"
params.bam = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1.pri.bam"
params.bam_index = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1.pri.bam.bai"
params.vcf_input = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1_deepvariant.vcf.gz"
params.vcf_index = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1_deepvariant.vcf.gz.tbi"

def get_sequence_map(fai_file) {
    def n_sequences = 0
    def max_length = 0
    def total_length = 0
    fai_file.eachLine { line ->
        def lspl   = line.split('\t')
        def chrom  = lspl[0]
        def length = lspl[1].toInteger()
        n_sequences ++
        total_length += length
        if (length > max_length) {
            max_length = length
        }
    }

    def sequence_map = [:]
    sequence_map.n_sequences = n_sequences
    sequence_map.total_length = total_length
    if (n_sequences) {
        sequence_map.max_length = max_length
    }
    return sequence_map
}

//
// Run Himut
//

include { HIMUT                    }   from '../../modules/local/himut.nf'
include { TABIX_BGZIP as BGZIP     }   from '../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX as TABIX_CSI }   from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_TBI }   from '../../modules/nf-core/tabix/tabix/main'

workflow {
    // RUN_HIMUT {

    // take:
    // ch_genome,              // genome fasta
    // ch_genome_index_fai,    // genome index
    // ch_region_list,         // region list, generated from .fai file
    // ch_aligned_reads,       // bam
    // ch_aligned_reads_index, // bam index
    // vcf_input,              // vcf input
    // vcf_index               // vcf index

    // main:
    // ch_versions = Channel.empty()


    // Run HIMUT

    // HIMUT (
    //     ch_genome,
    //     ch_region_list,
    //     ch_aligned_reads,
    //     ch_aligned_reads_index,
    //     vcf_input,
    //     vcf_index
    // )
    // ch_versions = ch_versions.mix ( HIMUT.out.versions.first() )

    HIMUT(params.fasta,
        params.fasta_index,
        params.assembly_report,
        params.bam,
        params.bam_index,
        params.vcf_input,
        params.vcf_index)

    ch_versions = Channel.empty()

    // compress the himut vcf files
    ch_meta = Channel.of([id:"test"])
    ch_meta
        .combine(HIMUT.out.vcf_output)
        .set{ch_himut_vcf}
    // ch_meta.view()

    ch_compressed_himut_vcf = BGZIP ( ch_himut_vcf ).output
    // ch_compressed_himut_vcf = BGZIP ( HIMUT.out.vcf_output ).output
    ch_versions = Channel.empty()
    // ch_versions             = ch_versions.mix ( BGZIP.out.versions.first() )


    ch_meta = Channel.of([id:"test"])
    ch_fasta_index = Channel.fromPath(params.fasta_index)
    // ch_fasta_index
    ch_meta
        .combine(ch_fasta_index)
        .map { meta, index -> [ [ id: meta.id ] + get_sequence_map(index) ] }
        .set { ch_genome_info }

    // index the compressed himut vcf files
    ch_compressed_himut_vcf
        .combine(ch_genome_info)
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
