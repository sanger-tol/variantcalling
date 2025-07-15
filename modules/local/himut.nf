params.fasta = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.fasta"
params.fasta_index = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.fasta.fai"
params.assembly_report = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1_assembly_report.txt"
// params.region_list = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/all_regions.GCA_937595015.1.txt"
params.bam = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1.pri.bam"
params.bam_index = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1.pri.bam.bai"
params.vcf_input = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1_deepvariant.vcf.gz"
// params.vcf_input = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1_deepvariant.vcf.bgz"
params.vcf_index = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1_deepvariant.vcf.gz.tbi"
// params.vcf_index = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets/data/GCA_937595015.1.pacbio.ilPolIcar1_deepvariant.vcf.bgz.tbi"

params.outdir = "/nfs/treeoflife-01/teams/tolit/users/yz12/pipelines/variant_calling/add_himut/assets"


process HIMUT {
    tag "$bam.baseName"
    // label 'process_single'
    // publishDir "${params.outdir}/himut", mode: 'copy', overwrite: true

    container "quay.io/sanger-tol/himut:1.0.0-c1"

    input:
    path fasta
    path fasta_index
    path assembly_report
    // path region_list
    path bam
    path bam_index
    path vcf_input //, stageAs: "input_vcf.bgz"
    path vcf_index //, stageAs: "input_vcf.bgz.tbi"
    // tuple val(meta), path(fasta), path(fasta_index)
    // tuple val(meta), path(bam), path(bam_index)
    // tuple val(meta), path(vcf_input), path(vcf_index)

    output:
    // tuple val(meta), path("*.vcf"), emit: vcf_output
    path  ("*.vcf"), emit: vcf_output
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    println vcf_input
    println vcf_index
    // --vcf ${vcf_input} \\

    """
    ln -s ${vcf_input} input.vcf.bgz
    ln -s ${vcf_index} input.vcf.bgz.tbi

    awk -F"\t" '\$3=="X" || \$3=="Y" || \$3=="Z" || \$3=="W" {print \$5}' ${assembly_report} > sex_chromosomes.txt
    cut -f1 ${fasta_index} | grep -vFxf sex_chromosomes.txt > regions_list.txt

    himut call \\
        -i ${bam} \\
        --ref ${fasta} \\
        --region_list regions_list.txt \\
        --vcf input.vcf.bgz \\
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
        params.fasta_index,
        params.assembly_report,
        // params.region_list,
        params.bam,
        params.bam_index,
        params.vcf_input,
        params.vcf_index)
}
