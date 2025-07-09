// // If region list is not provided, generate region_list from .fai file
// // Do not include sex chromosomes
// if (params.region_list) {
//     ch_region_list = Channel.fromPath(params.region_list)
// } else {
//     def make_region_list = { meta, fai ->
//         def chr_list = []
//         def chr_length = 0
//         def chr_name = ''
//         fai.splitEachLine('\t') { line ->
//             // How to select out sex chromosomes?
//             // if (line[0].startsWith('chr') && !line[0].contains('X') && !line[0].contains('Y') && !line[0].contains('Z') && !line[0].contains('W')) {
//                 chr_name = line[0]
//                 chr_length = line[1] as int
//                 if (chr_length > 0) {
//                     chr_list << chr_name
//                 }
//             //}
//         }
//         return [meta, chr_list]

//         ch_genome_index_fai
//             .map { meta, fai -> make_region_list(meta, fai) }
//             .set { ch_region_list }
//     }
// }

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
    if [param.region_list] ; then
        region_list=${region_list}
    else
        # Generate region list from .fai file
        region_list=\$(mktemp --suffix=.txt)
        grep -v -E '^(chrX|chrY|chrZ|chrW)' ${fasta}.fai | cut -f1 > \${region_list}
    fi

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