process HIMUT {
    // tag "$meta.id"
    tag "$bam.baseName"
    label 'process_medium'

    // container "quay.io/sanger-tol/himut:1.0.0-c1"
    // container "docker.io/sainsachiko/himut:1.0.0"
    container "docker.io/sainsachiko/himut:1.0.0-test"

    input:
    // tuple val(meta), path(fasta)
    // tuple val(meta), path(fasta_index)
    // tuple val(meta), path(assembly_report)
    // tuple val(meta), path(bam)
    // tuple val(meta), path(bam_index)
    // tuple val(meta), path(vcf_input)
    // tuple val(meta), path(vcf_index)

    // path(fasta)
    // path(fasta_index)
    // path(assembly_report)
    // path(bam)
    // path(bam_index)
    // path(vcf_input)
    // path(vcf_index)


    tuple val(meta), path(fasta)
    tuple val(meta), path(fasta_index)
    path(assembly_report)

    tuple val(meta), path(bam)
    // path(bam)
    tuple val(meta), path(bam_index)
    // path(bam_index)

    tuple val(meta), path(vcf_input)
    tuple val(meta), path(vcf_index)


    output:
    // tuple val(meta), path("*.vcf") , emit: vcf_output
    path  "*.vcf"         , emit: vcf_output
    path  "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Temporary rename the compressed VCF files (.gz) as input.vcf.bgz, as the vcf_input of Himut
    // Identify sex chromosomes (labelled with X/Y/Z/W in assembly report column 3), generate a sex_chromosomes.txt that contains the respective accession numbers
    // List the first column of fasta index (all the chromosome accession numbers), exclude those in sex_chromosomes.txt,
    //     generate region_list.txt as the region_list input of Himut
    """
    ln -s ${bam} input.bam
    ln -s ${bam_index} input.bam.bai
    ln -s ${vcf_input} input.vcf.bgz
    ln -s ${vcf_index} input.vcf.bgz.tbi

    awk -F"\t" '\$3=="X" || \$3=="Y" || \$3=="Z" || \$3=="W" {print \$5}' ${assembly_report} > sex_chromosomes.txt
    cut -f1 ${fasta_index} | grep -vFxf sex_chromosomes.txt > regions_list.txt

    himut call \\
        -i input.bam \\
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
