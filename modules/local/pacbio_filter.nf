process PACBIO_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gawk:5.1.0'
        : 'quay.io/biocontainers/gawk:5.1.0'}"

    input:
    tuple val(meta), path(txt)

    output:
    tuple val(meta), path("*.blocklist"), emit: list
    tuple val("${task.process}"), val('gawk'), eval("awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//'"), topic: versions, emit: versions_gawk

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pacbio_filter.sh ${txt} ${prefix}.blocklist
    """
}
