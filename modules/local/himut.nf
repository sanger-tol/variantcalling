process HIMUT {
    tag "$meta.id"
    // label 'process_high'

    container "quay.io/sanger-tol/himut:1.0.0-c2"

    input:

    tuple val(meta), path(fasta)
    tuple val(meta), path(fasta_index)
    path (assembly_report)
    tuple val(meta), path(bam), path(bam_index)
    tuple val(meta), path(vcf_input)
    tuple val(meta), path(vcf_index)


    output:
    tuple val(meta), path("*.somatic.vcf")                           , emit: vcf_output
    tuple val(meta), path("*.somatic.single_molecule_mutations.vcf") , emit: smm_vcf_output
    path  "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${bam.baseName}"

    // Himut takes reference fasta, fai, regions list, unaligned bam, bai, and germline vcf file as the input,
    //     and output somatic vcf (alternate count >= 1) + single_molecule_mutation somatic vcf (alternate count = 1).
    // Input requirements are:
    //     (1) Pair of fasta/fai and bam/bai files need to have the same name prefix
    //     (2) vcf file in .bgz format
    // So here, the codes temporary rename bam/bai and the vcf.gz/vcf_index, as the input of Himut.

    // To exclude sex chromosomes (labelled with X/Y/Z/W in assembly report column 3) from the analysis,
    //     first generate a sex_chromosomes.txt that contains the respective accession numbers (if no sex chromosomes identified, the txt file is empty);
    // List the first column of fai file (all the chromosome accession numbers), exclude those in sex_chromosomes.txt,
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
        -o ${prefix}.somatic.vcf \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        himut: \$(himut --version | sed 's/^himut //; s/ .*\$//')
    END_VERSIONS
    """
}
