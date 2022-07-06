////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////
nextflow.enable.dsl = 2

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, params.multiqc_config
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// channel to input file
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!'}

// Check segments
def segments  = []
def segmentsList = ['S', 'M', 'L']
if (!segmentsList.contains(params.segment)) {
    exit 1, "Invalid segment option: ${params.segment}. Valid options: ${segmentsList.join(', ')}"
} else {
    ch_prefix = params.segment + '-Segment'
    segments << params.segment
}

// Check primer scheme version
def primerschemes = []
def primerschemesList = ['V1', 'V2', 'V3']
if (!primerschemesList.contains(params.primer_scheme_version)) {
    exit 1, "Invalid primer scheme version option: ${params.primer_scheme_version}. Valid options: ${primerschemesList.join(', ')}"
} else {
    primerschemes << params.primer_scheme_version
}

// Check trimming parameters
def trimmersList = ['trimmomatic', 'fastp']
if (!params.skip_trimming) {
    if (!trimmersList.contains(params.trimmer)) {
        exit 1, "Invalid trimmer option: ${params.trimmer}. Valid options: ${trimmersList.join(', ')}"
    }
}

// Check alignment parameters
def prepareToolIndices  = []
def alignerList         = ['bwa', 'bowtie2']
if (!params.skip_alignment) {
    if (!alignerList.contains(params.aligner)) {
        exit 1, "Invalid aligner option: ${params.aligner}. Valid options: ${alignerList.join(', ')}"
    }
    prepareToolIndices << params.aligner
}

// check metadata file
if (!params.metadata) {
    ch_metadata = file("$projectDir/assets/metadata.csv", checkIfExists: true)
} else {
    ch_metadata = file(params.metadata)
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()



////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]

// add options to bwa index module
def bwa_index_options     = modules['bwa_index']
if (!params.save_reference)  { bwa_index_options['publish_files'] = false }

// add options to bwa mem module
def bwa_mem_options         = modules['bwa_mem']
if (params.save_align_intermeds) { bwa_mem_options.publish_files.put('bam','') }

// add options to bowtie2 build module
def bowtie2_build_options     = modules['bowtie2_build']
if (!params.save_reference)  { bowtie2_build_options['publish_files'] = false }

// add options to bowtie2 align module
def bowtie2_align_options         = modules['bowtie2_align']
if (params.save_align_intermeds) { bowtie2_align_options.publish_files.put('bam','') }

// iVar trim options
def ivar_trim_options            = modules['ivar_trim']
ivar_trim_options.args           += params.ivar_trim_noprimer ? "" : " -e"

// add options to samtools stats
def samtools_sort_options = modules['samtools_sort']
if (['bwa','bowtie2'].contains(params.aligner)) {
    if (params.save_align_intermeds && params.skip_markduplicates) {
        samtools_sort_options.publish_files.put('bam','')
        samtools_sort_options.publish_files.put('bai','')
    }
}
def samtools_sort_options_host = modules['samtools_sort_host']
if (['bwa','bowtie2'].contains(params.aligner)) {
    if (params.save_align_intermeds && params.skip_markduplicates) {
        samtools_sort_options.publish_files.put('bam','')
        samtools_sort_options.publish_files.put('bai','')
    }
}

def samtools_sort_options_reference = modules['samtools_sort_reference']
if (['bwa','bowtie2'].contains(params.aligner)) {
    if (params.save_align_intermeds && params.skip_markduplicates) {
        samtools_sort_options.publish_files.put('bam','')
        samtools_sort_options.publish_files.put('bai','')
    }
}

// Get total number of mapped reads from flagstat file
def get_mapped_from_flagstat(flagstat) {
    def mapped = 0
    flagstat.eachLine { line ->
        if (line.contains(' mapped (')) {
            mapped = line.tokenize().first().toInteger()
        }
    }
    return mapped
}


// Function that checks the number of mapped reads from flagstat output and returns true if > params.min_mapped_reads and otherwise false
pass_mapped_reads = [:]
fail_mapped_reads = [:]
def check_mapped(sample, flagstat, min_mapped=200) {
    mapped = get_mapped_from_flagstat(flagstat)
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    if (mapped < min_mapped.toInteger()) {
        log.info ">${c_red}>>>> $sample FAILED MAPPED READ THRESHOLD: ${mapped} < ${params.min_mapped}. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! <<<<${c_reset}<"
        fail_mapped_reads[sample] = mapped
        return false
    } else {
        pass_mapped_reads[sample] = mapped
        return true
    }
}


def count_lines(vcf) {
    def lines = vcf.countLines()
    return lines
}

pass_vcfs = [:]
fail_vcfs = [:]
def check_vcf(sample, vcf, header=14) {
    lines = count_lines(vcf)
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    if (lines == header.toInteger()) {
        log.info ">${c_red}>>>> $sample FAILED VCF THRESHOLD: ${lines} == ${header}. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! <<<<${c_reset}<"
        fail_vcfs[sample] = lines
        return false
    } else {
        pass_vcfs[sample] = lines
        return true
    }
}

// add options to multiqc module
def multiqc_options         = modules['multiqc']
multiqc_options.args       += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

include { readInputFile                               } from './modules/local/process/samplesheet_check'                     addParams( options: [:]                                                                                                                                                                        )
include { MULTIQC                                     } from './modules/local/process/multiqc'                               addParams( options: multiqc_options                                                                                                                                                            )
include { GET_SOFTWARE_VERSIONS                       } from './modules/local/process/get_software_versions'                 addParams( options: [publish_files : ['csv':'']]                                                                                                                                               )
include { TRIMMOMATIC                                 } from './modules/local/process/trimmomatic'                           addParams( options: [:]                                                                                                                                                                        )
include { PREPARE_GENOME                              } from './modules/local/subworkflow/prepare_genome'                    addParams( genome_options: publish_genome_options, index_options: publish_index_options, bwa_index_options: bwa_index_options, bowtie2_index_options: bowtie2_build_options                    )
include { PREPARE_PRIMER_SCHEMES                      } from './modules/local/subworkflow/prepare_primer_schemes'            addParams( options: [:]                                                                                                                                                                        )
include { PREPARE_HOST_GENOME                         } from './modules/local/subworkflow/prepare_host_genome'               addParams( genome_options: publish_genome_options, index_options: publish_index_options, bwa_index_options: bwa_index_options, bowtie2_index_options: bowtie2_build_options                    )
include { ALIGN_BWA as ALIGN_BWA_HOST                 } from './modules/local/subworkflow/align_bwa'                         addParams( bwa_mem_options: bwa_mem_options, samtools_options: samtools_sort_options_host                                                                                                      )
include { ALIGN_BWA as ALIGN_BWA_FASTA                } from './modules/local/subworkflow/align_bwa'                         addParams( align_options:   bwa_mem_options, samtools_options: samtools_sort_options_reference                                                                                                 )
include { ALIGN_BOWTIE2 as ALIGN_BOWTIE2_HOST         } from './modules/local/subworkflow/align_bowtie2'                     addParams( align_options: bowtie2_align_options, samtools_options: samtools_sort_options_host                                                                                                  )
include { ALIGN_BOWTIE2 as ALIGN_BOWTIE2_FASTA        } from './modules/local/subworkflow/align_bowtie2'                     addParams( align_options: bowtie2_align_options, samtools_options: samtools_sort_options_reference                                                                                             )
include { SAMTOOLS_UNMAPPED_FASTQ as HOST_UNMAPPED    } from './modules/local/process/samtools_unmapped_fastq'               addParams( options: modules['samtools_unmapped_fastq_host']                                                                                                                                    )
include { PARSE_FLAGSTAT as PARSE_ALIGNMENT_FLAGSTAT  } from './modules/local/process/parse_alignment_stats'                 addParams( options: modules['alignment_stats']                                                                                                                                                 )
include { PARSE_FLAGSTAT as PARSE_IVAR_FLAGSTAT       } from './modules/local/process/parse_alignment_stats'                 addParams( options: modules['ivar_trim']                                                                                                                                                       )
include { SAMTOOLS_MAPPED_BAM                         } from './modules/local/process/samtools_mapped_bam'                   addParams( options: modules['samtools_mapped_bam']                                                                                                                                             )
include { MARKDUPLICATES_PICARD                       } from './modules/local/subworkflow/markduplicates_picard'             addParams( picard_markduplicates_options: modules['picard_markduplicates'], samtools_options: modules['picard_markduplicates_samtools']                                                        )
include { TRIM_IVAR                                   } from './modules/local/subworkflow/trim_ivar'                         addParams( options: modules['ivar_trim'], ivar_trim_options: ivar_trim_options, samtools_options: modules['ivar_trim_samtools']                                                                )
include { SNPEFF_ANNOTATE                             } from './modules/local/process/snpeff_annotate'                       addParams(options: [:]                                                                                                                                                                         )
include { BEDTOOLS_GENOME_COVERAGE                    } from './modules/local/process/bedtools_genome_coverage'              addParams( options: modules['bedtools_genomecov']                                                                                                                                              )
include { BEDTOOLS_AMPLICON_COVERAGE                  } from './modules/local/process/bedtools_amplicon_coverage'            addParams( options: modules['bedtools_ampliconcov']                                                                                                                                            )
include { BEDTOOLS_AMPLICON_MEAN_COVERAGE             } from './modules/local/process/bedtools_amplicon_mean_coverage'       addParams( options: modules['bedtools_amplicon_mean_coverage']                                                                                                                                 )
include { BEDTOOLS_AMPLICON_BASE_COVERAGE             } from './modules/local/process/bedtools_amplicon_base_coverage'       addParams( options: modules['bedtools_amplicon_base_coverage']                                                                                                                                 )
include { COVERAGE_PLOTS                              } from './modules/local/process/coverage_plots'                        addParams( options: modules['plots']                                                                                                                                                           )
include { CONSENSUS_QC                                } from './modules/local/process/consensus_qc'                          addParams( options: modules['consensus_qc']                                                                                                                                                    )
include { QC_SUMMARY_CSV                              } from './modules/local/process/qc_summary_csv'                        addParams( options: modules['consensus_qc']                                                                                                                                                    )
include { READ_COUNT as READS_RAW                     } from './modules/local/process/read_count'                            addParams( options: modules['raw_counts']                                                                                                                                                      )
include { READ_COUNT as READS_TRIMMED                 } from './modules/local/process/read_count'                            addParams( options: modules['trimmed_counts']                                                                                                                                                  )
include { MERGE_COUNTS as MERGE_RAW_COUNTS            } from './modules/local/process/merge_counts'                          addParams( options: modules['merge_raw_counts']                                                                                                                                                )
include { MERGE_COUNTS as MERGE_TRIMMED_COUNTS        } from './modules/local/process/merge_counts'                          addParams( options: modules['merge_trimmed_counts']                                                                                                                                            )
include { MERGE_COUNTS as MERGE_MAPPED_COUNTS         } from './modules/local/process/merge_counts'                          addParams( options: modules['merge_mapped_counts']                                                                                                                                             )
include { MERGE_COUNTS as MERGE_PRIMERTRIMMED         } from './modules/local/process/merge_counts'                          addParams( options: modules['merge_primmermapped']                                                                                                                                             )

include { MERGE_STATS                                 } from './modules/local/process/merge_stats'                           addParams( options: modules['summary_stats']                                                                                                                                                   )
include { PLOT_READ_COUNTS                            } from './modules/local/process/plot_read_counts'                      addParams( options: modules['plots']                                                                                                                                                           )


////////////////////////////////////////////////////
/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
////////////////////////////////////////////////////

/*
 * MODULE: Installed directly from nf-core/modules
 */

include { FASTQC                     } from './modules/nf-core/software/fastqc/main'                 addParams( options: modules['fastqc']          )
include { FASTP                      } from './modules/nf-core/software/fastp/main'                  addParams( options: modules['fastp']           )
include { IVAR_VARIANTS              } from './modules/nf-core/software/ivar/variants/main'          addParams( options: modules['ivar_variants']   )
include { IVAR_CONSENSUS             } from './modules/nf-core/software/ivar/consensus/main'         addParams( options: modules['ivar_consensus']  )


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []

workflow RVFVAMPLICONSEQ {

    ch_software_versions = Channel.empty()

    // SUBWORKFLOW: Read samplesheet and stage input reads
    ch_reads = Channel.empty()
    ch_reads = readInputFile(ch_input, params.single_end)
    
    // MODULE: FastQC
    ch_fastqc_multiqc = Channel.empty()
    if (!params.skip_qc) {
        FASTQC (ch_reads)
        ch_fastqc_multiqc    = FASTQC.out.zip
        ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

        // count raw reads
        READS_RAW(ch_reads)
        ch_raw_counts = READS_RAW.out.csv
    }

    // MODULE: Trim reads
    if (params.skip_trimming) {
        ch_trimmed_reads     = ch_reads
        ch_trimmed_counts    = ch_raw_counts
    }

    ch_fastp_multiqc = Channel.empty()
    if (!params.skip_trimming && params.trimmer == 'fastp') {
        FASTP (ch_reads)
        ch_trimmed_reads     = FASTP.out.reads
        ch_fastp_multiqc     = FASTP.out.log
        ch_software_versions = ch_software_versions.mix(FASTP.out.version.first().ifEmpty(null))

        // count trimmed reads
        READS_TRIMMED(ch_trimmed_reads)
        ch_trimmed_counts    = READS_TRIMMED.out.csv
    }
    
    ch_trimmomatic_multiqc = Channel.empty()
    if (!params.skip_trimming && params.trimmer == 'trimmomatic' ) {
        TRIMMOMATIC(ch_reads)
        ch_trimmed_reads       = TRIMMOMATIC.out.reads
        ch_trimmomatic_multiqc = TRIMMOMATIC.out.log
        ch_software_versions   = ch_software_versions.mix(TRIMMOMATIC.out.version.first().ifEmpty(null))

        // count trimmed reads
        ch_trimmed_counts    = READS_TRIMMED(ch_trimmed_reads)
    }

    // SUBWORKFLOW: Filter host reads and extract unmapped reads
    if (!params.host_fasta || !params.host_bwa_index || !params.host_bowtie2_index) {
        ch_bwa_multiqc          = Channel.empty()
        ch_bowtie2_multiqc      = Channel.empty()
        ch_host_filtered_reads  = ch_trimmed_reads
    }
    if (params.host_fasta && params.aligner == 'bwa') {
       PREPARE_HOST_GENOME(prepareToolIndices)
        ALIGN_BWA_HOST (
            ch_trimmed_reads,
            PREPARE_HOST_GENOME.out.bwa_index
        )
        ch_genome_bam           = ALIGN_BWA_HOST.out.bam
        ch_genome_bai           = ALIGN_BWA_HOST.out.bai
        ch_samtools_stats       = ALIGN_BWA_HOST.out.stats
        ch_samtools_flagstat    = ALIGN_BWA_HOST.out.flagstat
        ch_samtools_idxstats    = ALIGN_BWA_HOST.out.idxstats
        HOST_UNMAPPED(ch_genome_bam)
        ch_host_filtered_reads  = HOST_UNMAPPED.out.reads
    }
    if (params.host_fasta && params.aligner == 'bowtie2') {
        PREPARE_HOST_GENOME(prepareToolIndices)
        ALIGN_BOWTIE2_HOST (
            ch_trimmed_reads,
            PREPARE_HOST_GENOME.out.bowtie2_index.collect()
        )
        ch_genome_bam        = ALIGN_BOWTIE2_HOST.out.bam
        ch_genome_bai        = ALIGN_BOWTIE2_HOST.out.bai
        ch_bowtie2_multiqc   = ALIGN_BOWTIE2_HOST.out.log
        ch_samtools_stats    = ALIGN_BOWTIE2_HOST.out.stats
        ch_samtools_flagstat = ALIGN_BOWTIE2_HOST.out.flagstat
        ch_samtools_idxstats = ALIGN_BOWTIE2_HOST.out.idxstats
        HOST_UNMAPPED(ch_genome_bam)
        ch_host_filtered_reads  = HOST_UNMAPPED.out.reads
    }
    if (params.host_bwa_index) {
        PREPARE_HOST_GENOME(prepareToolIndices)
        ALIGN_BWA_HOST (
            ch_trimmed_reads,
            PREPARE_HOST_GENOME.out.bwa_index
        )
        ch_genome_bam           = ALIGN_BWA_HOST.out.bam
        ch_genome_bai           = ALIGN_BWA_HOST.out.bai
        ch_samtools_stats       = ALIGN_BWA_HOST.out.stats
        ch_samtools_flagstat    = ALIGN_BWA_HOST.out.flagstat
        ch_samtools_idxstats    = ALIGN_BWA_HOST.out.idxstats
        HOST_UNMAPPED(ch_genome_bam)
        ch_host_filtered_reads  = HOST_UNMAPPED.out.reads
    } 
    if (params.host_bowtie2_index) {
        PREPARE_HOST_GENOME(prepareToolIndices)
        ALIGN_BOWTIE2_HOST (
            ch_trimmed_reads,
            PREPARE_HOST_GENOME.out.bowtie2_index
        )
        ch_genome_bam           = ALIGN_BOWTIE2_HOST.out.bam
        ch_genome_bai           = ALIGN_BOWTIE2_HOST.out.bai
        ch_bowtie2_multiqc      = ALIGN_BOWTIE2_HOST.out.log
        ch_samtools_stats       = ALIGN_BOWTIE2_HOST.out.stats
        ch_samtools_flagstat    = ALIGN_BOWTIE2_HOST.out.flagstat
        ch_samtools_idxstats    = ALIGN_BOWTIE2_HOST.out.idxstats
        HOST_UNMAPPED(ch_genome_bam)
        ch_host_filtered_reads  = HOST_UNMAPPED.out.reads
    }


    // SUBWORKFLOW: Alignment to the reference FASTA genome with BWA MEM
    ch_bwa_multiqc = Channel.empty()
    ch_sort_bam    = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'bwa') {        
        PREPARE_GENOME(prepareToolIndices)
        ALIGN_BWA_FASTA (
            ch_host_filtered_reads,
            PREPARE_GENOME.out.bwa_index
        )
        ch_fasta             = PREPARE_GENOME.out.fasta
        ch_gff               = PREPARE_GENOME.out.gff
        ch_genome_bam        = ALIGN_BWA_FASTA.out.bam
        ch_genome_bai        = ALIGN_BWA_FASTA.out.bai
        ch_samtools_stats    = ALIGN_BWA_FASTA.out.stats
        ch_samtools_flagstat = ALIGN_BWA_FASTA.out.flagstat
        ch_samtools_idxstats = ALIGN_BWA_FASTA.out.idxstats
    }

    // SUBWORKFLOW: Alignment to the reference FASTA genome with BOWTIE2
    ch_bowtie2_multiqc = Channel.empty()
    if (!params.skip_alignment && params.aligner == 'bowtie2') {    
        PREPARE_GENOME(prepareToolIndices)    
        ALIGN_BOWTIE2_FASTA (
            ch_host_filtered_reads,
            PREPARE_GENOME.out.bowtie2_index
        )
        ch_fasta             = PREPARE_GENOME.out.fasta
        ch_gff               = PREPARE_GENOME.out.gff
        ch_genome_bam        = ALIGN_BOWTIE2_FASTA.out.bam
        ch_genome_bai        = ALIGN_BOWTIE2_FASTA.out.bai
        ch_bowtie2_multiqc   = ALIGN_BOWTIE2_FASTA.out.log
        ch_samtools_stats    = ALIGN_BOWTIE2_FASTA.out.stats
        ch_samtools_flagstat = ALIGN_BOWTIE2_FASTA.out.flagstat
        ch_samtools_idxstats = ALIGN_BOWTIE2_FASTA.out.idxstats
    }

    // add flagstat channel
    ch_sort_bam = ch_genome_bam.join(ch_samtools_flagstat)
    ch_alignment_stats = PARSE_ALIGNMENT_FLAGSTAT(ch_samtools_flagstat)

    // MODULE: merge the read counts and mapping stats
    
    ch_raw_counts
        .filter { row -> file(row[1]) }
        .map { it[1] }
        .set { ch_raw_counts }

    ch_trimmed_counts
        .filter { row -> file(row[1]) }
        .map { it[1] }
        .set { ch_trimmed_counts }

    ch_alignment_stats
        .filter { row -> file(row[1]) }
        .map { it[1] }
        .set { ch_alignment_stats }

    MERGE_RAW_COUNTS ( ch_raw_counts.collect(), ch_prefix )
    MERGE_TRIMMED_COUNTS ( ch_trimmed_counts.collect(), ch_prefix )
    MERGE_MAPPED_COUNTS ( ch_alignment_stats.collect(), ch_prefix )

    ch_raw_rc = MERGE_RAW_COUNTS.out.csv
    ch_trimmed_rc = MERGE_TRIMMED_COUNTS.out.csv
    ch_mapped_rc  = MERGE_MAPPED_COUNTS.out.csv

    MERGE_STATS (ch_raw_rc, ch_trimmed_rc, ch_mapped_rc, ch_prefix)
    PLOT_READ_COUNTS (MERGE_STATS.out.csv, ch_prefix)

    // Remove samples that failed mapped read threshold
    ch_sort_bam
        .filter { row -> check_mapped( row[0], row[2], params.min_mapped )  }
        .map {it[0..1]}
        .set { ch_sort_bam }

    // SUBWORKFLOW: Extract mapped alignments and index
    SAMTOOLS_MAPPED_BAM (ch_sort_bam)


    // SUBWORKFLOW: prepare primer schemes
    PREPARE_PRIMER_SCHEMES(segments, primerschemes)

    // SUBWORKFLOW: Trim amplicon primers 
    TRIM_IVAR (
        SAMTOOLS_MAPPED_BAM.out.bam.join(SAMTOOLS_MAPPED_BAM.out.bai, by: [0]),
        PREPARE_PRIMER_SCHEMES.out.bed.collect()
    )
    ch_ivar_trim_stats = PARSE_IVAR_FLAGSTAT(TRIM_IVAR.out.flagstat)

    ch_ivar_trim_stats
        .filter { row -> file(row[1]) }
        .map { it[1] }
        .set { ch_ivar_trim_stats }


    // MODULE: merge alignment metrics and summarize with read counts

    // MERGE_PRIMERTRIMMED ( ch_ivar_trim_stats.collect(), ch_prefix )
    // ch_primertrimmed = MERGE_PRIMERTRIMMED.out.csv
    // MERGE_STATS (ch_raw_rc, ch_trimmed_rc, ch_mapped_rc, ch_primertrimmed, ch_prefix)

    // SUBWORKFLOW: Mark duplicates
    ch_ivar_trim_bam = TRIM_IVAR.out.bam
    if (!params.skip_markduplicates) {
        MARKDUPLICATES_PICARD(ch_ivar_trim_bam)
        ch_markduplicates_bam = MARKDUPLICATES_PICARD.out.bam
    } else {
        ch_markduplicates_bam = ch_ivar_trim_bam
    }

    // MODULE: Genome wide coverage
    BEDTOOLS_GENOME_COVERAGE (ch_markduplicates_bam)

    //MODULE: AMPLICON coverage
    BEDTOOLS_AMPLICON_COVERAGE (ch_markduplicates_bam, PREPARE_PRIMER_SCHEMES.out.bed.collect())

    // MODULE: AMPLICON mean coverage
    BEDTOOLS_AMPLICON_MEAN_COVERAGE (ch_markduplicates_bam, PREPARE_PRIMER_SCHEMES.out.bed.collect())

    // MODULE: AMPLICON base coverage
    BEDTOOLS_AMPLICON_BASE_COVERAGE (ch_markduplicates_bam, PREPARE_PRIMER_SCHEMES.out.bed.collect())

    // MODULE: call variants with iVar
    IVAR_VARIANTS (ch_markduplicates_bam, ch_fasta.collect(), ch_gff.collect())

    // MODULE: generate consensus sequence with iVar
    IVAR_CONSENSUS (ch_markduplicates_bam, ch_fasta.collect())

    // MODULE: generate consensus qc metrics
    ch_consensus_qc = ch_markduplicates_bam.join(IVAR_CONSENSUS.out.fasta)
    CONSENSUS_QC (ch_consensus_qc, ch_fasta.collect())



    // Remove samples that failed vcf
    ch_vcf = IVAR_VARIANTS.out.vcf

    ch_vcf
        .filter { row -> check_vcf( row[0], row[1], 14 )  }
        .map {it[0..1]}
        .set { ch_vcf }

    // MODULE: Build SnpEff db and Annotate variants
    SNPEFF_ANNOTATE(ch_vcf, ch_fasta.collect(), ch_gff.collect())

    // MODULE: merge qc csv files
    ch_qc_csv = CONSENSUS_QC.out.csv
    ch_qc_csv
        .filter { row -> file(row[1])}
        .map { it[1] }
        .set { ch_qc_csv }

    QC_SUMMARY_CSV (ch_qc_csv.collect(), MERGE_STATS.out.csv, ch_prefix)

    // channel for all the consensus fasta files
    ch_ivar_consensus_fasta = IVAR_CONSENSUS.out.fasta
    ch_ivar_consensus_fasta
        .filter { row -> file(row[1])}
        .map { it[1] }
        .set { ch_ivar_consensus_fasta }

    // channel for all the genomewide coverage files
    ch_genome_coverage_bed = BEDTOOLS_GENOME_COVERAGE.out.coverage_bed
    ch_genome_coverage_bed
        .filter { row -> file(row[1])}
        .map { it[1] }
        .set { ch_genome_coverage_bed }


    // channel for all the amplicon per base coverage files
    ch_amplicon_mean_bed = BEDTOOLS_AMPLICON_MEAN_COVERAGE.out.amplicon_mean_coverage_bed
    ch_amplicon_mean_bed
        .filter { row -> file(row[1])}
        .map { it[1] }
        .set { ch_amplicon_mean_bed }


    // MODULE: Coverage plots
    COVERAGE_PLOTS (
        ch_ivar_consensus_fasta.collect(),
        ch_genome_coverage_bed.collect(),
        ch_amplicon_mean_bed.collect(),
        ch_metadata,
        QC_SUMMARY_CSV.out.csv,
        ch_prefix
    )


    // MODULE: Pipeline reporting
    
    GET_SOFTWARE_VERSIONS ( ch_software_versions.map { it }.collect())

    // MultiQC
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_fastqc_multiqc.collect{it[1]}.ifEmpty([]),
            ch_fastp_multiqc.collect{it[1]}.ifEmpty([]),
            //ch_bowtie2_multiqc.collect{it[1]}.ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////