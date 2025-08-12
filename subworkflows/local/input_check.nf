//
// Check input samplesheet and get read channels
//

workflow INPUT_CHECK {
    take:
    ch_samplesheet    // channel: [ val(meta), /path/to/reads ]

    main:
    ch_versions = Channel.empty()

    // ch_samplesheet.view()
    // Create the samplesheet channel
    ch_samplesheet
        .map { meta, datafile -> create_data_channel( meta, datafile ) }
        .set { reads }


    emit:
    reads                          // channel: [ val(meta), /path/to/datafile ]
    versions = ch_versions         // channel: [ versions.yml                 ]
}

// Function to get list of [ meta, VCF ]
def create_data_channel ( LinkedHashMap row, datafile ) {
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
