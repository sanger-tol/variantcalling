//
// Check input samplesheet and get a cram channel
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_cram_channel(it) }
        .set { cram }

    emit:
    cram                                     // channel: [ val(meta), [ cram ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ cram ] ]
def create_cram_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.datatype = row.datatype

    // add path(s) of the cram file to the meta map
    def cram_meta = []
    if (!file(row.datafile).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 cram file does not exist!\n${row.datafile}"
    }
    
    return cram_meta
}
