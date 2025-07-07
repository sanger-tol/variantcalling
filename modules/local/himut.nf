process HIMUT {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/repository/sanger-tol/himut"   // HIMUT version 1.0.0

    input:
    tuple val(meta), path(fasta)
    tuple val(meta), path(region_list)
    tuple val(meta), path(bam), path(bam_index)
    tuple val(meta), path(vcf_input), path(vcf_index)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf_output
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    himut call \\
        -i=${bam} \\
        --ref=${fasta} \\
        --region_list=${region_list} \\
        --vcf_input=${vcf_input} \\
        --vcf_index=${vcf_index} \\
        -o=${meta.id}.vcf \\
        --non_human_sample \\
        -t=${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        himut: \$(himut --version | sed 's/^himut //; s/ .*\$//')
    END_VERSIONS
    """
}