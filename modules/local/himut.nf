params.fasta = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/yz12-add_himut/assets/data/GCA_937595015.1.fasta"
// params.fasta_index = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/yz12-add_himut/assets/data/GCA_937595015.1.fasta.fai"
params.region_list = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/yz12-add_himut/assets/data/all_regions.GCA_937595015.1.txt"
params.bam = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/yz12-add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1.pri.bam"
params.bam_index = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/yz12-add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1.pri.bam.bai"
params.vcf_input = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/yz12-add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1_deepvariant.vcf.gz"
params.vcf_index = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/yz12-add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1_deepvariant.vcf.bgz.tbi"

params.outdir = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/yz12-add_himut/assets"


process HIMUT {
    tag "$bam.baseName"
    // label 'process_single'
    // publishDir "${params.outdir}/himut", mode: 'copy', overwrite: true

    container "quay.io/sanger-tol/himut:1.0.0-c1"

    input:
    path fasta
    // path fasta_index
    path region_list
    path bam
    path bam_index
    path vcf_input, stageAs( {bam.BaseName}.vcf.bgz )
    path vcf_index
    // tuple val(meta), path(fasta), path(fasta_index)
    // tuple val(meta), path(bam), path(bam_index)
    // tuple val(meta), path(vcf_input), path(vcf_index)

    output:
    // tuple path("*.vcf"), emit: vcf_output
    path  ("*.vcf"), emit: vcf_output
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    // if (params.region_list) {
    // ch_region_list = Channel.fromPath(params.region_list)
    // } else {
    // def make_region_list = { meta, fai ->
    //     def chr_list = []
    //     def chr_length = 0
    //     def chr_name = ''
    //     fai.splitEachLine('\t') { line ->
    //         // How to select out sex chromosomes?
    //         // if (line[0].startsWith('chr') && !line[0].contains('X') && !line[0].contains('Y') && !line[0].contains('Z') && !line[0].contains('W')) {
    //             chr_name = line[0]
    //             chr_length = line[1] as int
    //             if (chr_length > 0) {
    //                 chr_list << chr_name
    //             }
    //         //}
    //     }
    //     return [meta, chr_list]

    //     ch_genome_index_fai
    //         .map { meta, fai -> make_region_list(meta, fai) }
    //         .set { ch_region_list }
    // } }

    """
    himut call \\
        -i ${bam} \\
        --ref ${fasta} \\
        --region_list ${region_list} \\
        --vcf ${vcf_input} \\
        --non_human_sample \\
        -o ${bam.baseName}.himut.vcf \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        himut: \$(himut --version | sed 's/^himut //; s/ .*\$//')
    END_VERSIONS
    """
}

workflow {
    HIMUT(params.fasta,
        // params.fasta_index,
        params.region_list,
        params.bam,
        params.bam_index,
        params.vcf_input,
        params.vcf_index)
}
