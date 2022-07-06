#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/rvfvampliconseq
========================================================================================
 nf-core/rvfvampliconseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/rvfvampliconseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/rvfvampliconseq --input '*_R{1,2}.fastq.gz' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --        GENOME PARAMETER VALUES           -- */
////////////////////////////////////////////////////

// params.fasta         = Checks.get_genome_attribute(params, 'fasta')
// params.bed           = Checks.get_genome_attribute(params, 'bed')
// params.gtf           = Checks.get_genome_attribute(params, 'gtf')
// params.gff           = Checks.get_genome_attribute(params, 'gff')
// params.bwa_index     = Checks.get_genome_attribute(params, 'bwa')
// params.bowtie2_index = Checks.get_genome_attribute(params, 'bowtie2')


////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

// Check genome key exists if provided
Checks.genome_exists(params, log)


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
    include { RVFVAMPLICONSEQ } from './rvfvampliconseq' addParams( summary_params: summary_params )
    RVFVAMPLICONSEQ ()
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////


// if (params.validate_params) {
//     NfcoreSchema.validateParameters(params, json_schema, log)
// }

// ////////////////////////////////////////////////////
// /* --     Collect configuration parameters     -- */
// ////////////////////////////////////////////////////

// // Check if genome exists in the config file
// if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
//     exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(', ')}"
// }

// // TODO nf-core: Add any reference files that are needed
// // Configurable reference genomes
// //
// // NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// // If you want to use the channel below in a process, define the following:
// //   input:
// //   file fasta from ch_fasta
// //
// params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
// if (params.fasta) { ch_fasta = file(params.fasta, checkIfExists: true) }

// // Check AWS batch settings
// if (workflow.profile.contains('awsbatch')) {
//     // AWSBatch sanity checking
//     if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
//     // Check outdir paths to be S3 buckets if running on AWSBatch
//     // related: https://github.com/nextflow-io/nextflow/issues/813
//     if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
//     // Prevent trace files to be stored on S3 since S3 does not support rolling files.
//     if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
// }

// // Stage config files
// ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
// ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
// ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
// ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)