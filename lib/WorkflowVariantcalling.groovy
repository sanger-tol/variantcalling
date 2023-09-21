//
// This file holds several functions specific to the workflow/variantcalling.nf in the sanger-tol/variantcalling pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowVariantcalling {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        if (!params.fasta) {
            Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        }
    }
}
