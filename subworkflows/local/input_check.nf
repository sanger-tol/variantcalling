//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { [
            [ id: it.sample, sample: it.sample.replaceAll(/_T\d+$/, ''), type: it.datatype ], 
            file(it.datafile)
            ] }
        .set { reads }
        
    emit:
    reads                                     // channel: [ val(meta), data ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}
