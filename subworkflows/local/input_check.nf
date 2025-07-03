//
// Check input samplesheet and get read channels
//

include { samplesheetToList } from 'plugin/nf-schema'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    ch_versions = Channel.empty()

    Channel
        .fromList(samplesheetToList(samplesheet, "${projectDir}/assets/schema_input.json"))
        .map { check_data_channel( it ) }
        .set { reads }

    emit:
    reads                  // channel: [ val(meta), data ]
    versions = ch_versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, reads ]
def check_data_channel(meta, reads) {
    // create meta map
    def meta    = [:]
    meta.id     = row.sample
    meta.sample = row.sample.split('_')[0..-2].join('_')
    meta.datatype   = row.datatype

    if ( meta.datatype == "pacbio" ) {
        platform = "PACBIO"
    }
    meta.read_group  = "\'@RG\\tID:" + row.datafile.split('/')[-1].split('\\.')[0..-2].join('.') + "\\tPL:" + platform + "\\tSM:" + meta.sample + "\'"

    // add path(s) of the read file(s) to the meta map
    def data_meta = []
    if ( !file(row.datafile).exists() ) {
        exit 1, "ERROR: Please check input samplesheet -> Data file does not exist!\n${row.datafile}"
    } else {
        data_meta = [ meta, file(row.datafile) ]
    }
    return data_meta
}
