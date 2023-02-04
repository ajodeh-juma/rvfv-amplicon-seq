<!-- # ![nf-core/rvfvampliconseq](docs/images/nf-core-rvfvampliconseq_logo.png) -->

<!-- **Identify Rift Valley Fever virus lineages of a nucleotide sequence**. -->

  - [Introduction](#introduction)
  - [Installation](#installation)
  - [Testing](#test)
  - [Usage](#usage)
  - [Input](#input)
  - [Metadata](#metadata)
  - [Method details](#method-details)
  - [Output](#output)
  - [Pipeline summary](#pipeline-summary)
  - [Citations](#citations)
  

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23rvfvampliconseq-4A154B?logo=slack)](https://nfcore.slack.com/channels/rvfvampliconseq)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/ajodeh-juma/rvfvampliconseq/blob/master/LICENSE)
<!-- [![Docker](https://img.shields.io/docker/automated/nfcore/rvfvampliconseq.svg)](https://hub.docker.com/r/nfcore/rvfvampliconseq) -->
<!-- [![GitHub Actions CI Status](https://github.com/nf-core/rvfvampliconseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/rvfvampliconseq/actions) -->
<!-- [![GitHub Actions Linting Status](https://github.com/nf-core/rvfvampliconseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/rvfvampliconseq/actions) -->
<!-- ![](https://img.shields.io/badge/uses-docker-blue.svg) -->

[![Twitter Follow](https://img.shields.io/twitter/follow/john_juma.svg?style=social)](https://twitter.com/john_juma)

## Introduction

**rvfvampliconseq** A nextflow pipeline for analyzing Rift Valley fever virus amplicon sequencing data from Illumina instrument. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.


## Installation

**rvfvampliconseq** runs on UNIX/LINUX systems. You will install Miniconda3 from [here](https://docs.conda.io/en/latest/miniconda.html). Once Miniconda3 has been installed, proceed with pipeline installation

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline

```
git clone https://github.com/ajodeh-juma/rvfv-amplicon-seq.git
cd rvfv-amplicon-seq
conda env create -f environment.yml
conda activate rvfvampliconseq-env
```

## Testing

  - Optional: Test the installation on a minimal dataset bundled in the installation
    
    ```nextflow run main.nf -profile test```
    

## Usage

For minimal pipeline options, use the ```--help``` flag e.g. 

```nextflow run main.nf --help```

To see all the options, use the ```--show_hidden_params``` flag e.g.

```nextflow run main.nf --help --show_hidden_params```


## Input

The input is a comma-separated values (`CSV`) file having the columns: `sample`, `fastq_1` and `fastq_2` for paired-end sequence
data. The 'fastq_1' and 'fastq_2' columns should point to the absolute paths of the fastq files. The 'sample' name 
should correspond to the basename of the fastq files.

|__sample__ |__fastq_1__ |__fastq_2__ |
| --- | --- | --- |
|`AM-M1`|`/home/jjuma/data/genomics/rvfv/illumina/miseq/run/09032022/testdata/AM-M1_S1_L001_R1_001.fastq.gz`| `/home/jjuma/data/genomics/rvfv/illumina/miseq/run/09032022/testdata/AM-M1_S1_L001_R2_001.fastq.gz` |
|`RU-1`|`/Users/jjuma/data/genomics/rvfv/illumina/miseq/run/09032022/testdata/RU-1_S1_L001_R1_001.fastq.gz`| `/home/jjuma/data/genomics/rvfv/illumina/miseq/run/09032022/testdata/RU-1_S1_L001_R2_001.fastq.gz` |


## Metadata

If you intend to generate plots on coverage vs Ct values, include a metadata table in csv format having the columns
`sample_name`, `Ct`. For example:


|__sample_name__|__sample_type__|__host__|__platform__|__instrument__|__strategy__|__date__|__Ct__|__location__|__country__|culture|
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `AM-M1`|Aborted Foetus|Cow|Illumina|MiSeq|amplicon|16-01-2021|21.918|Kiambu|Kenya|No|
| `RU-1`|Serum|Cow|Illumina|MiSeq|amplicon|01-01-2018|25.245|Rulindo|Rwanda|No|


A typical command for S segment is
```
nextflow run main.nf \
      --input /home/jjuma/PhD_RVF2019/projects/RVFv_amplicon_tiling_PCR/qPCR_Read_outs/09032022_samplesheet.csv \
      --segment S \
      --skip_markduplicates \
      --metadata /home/jjuma/PhD_RVF2019/projects/RVFv_amplicon_tiling_PCR/qPCR_Read_outs/metadata_09032022.csv \
      --outdir "${OUTDIR}/S-segment-outdir" \
      -work-dir "${OUTDIR}/S-segment-workdir" \
      -resume
```
For the other segments, just replace the argument value for `--segment` to either `L` or `M`.
## Method details

The pipeline offers several parameters including as highlighted:

```
Input/output options
  --input                      [string]  Path to comma-separated file containing information about the samples in the experiment.
  --single_end                 [boolean] Specifies that the input is single-end reads.
  --outdir                     [string]  The output directory where the results will be saved. [default: ./results]
  --multiqc_title              [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.
  --email                      [string]  Email address for completion summary.

Reference genome options
  --segment                    [string]  genomic segment of the virus. options are 'S', 'M' and 'L'
  --host_fasta                 [string]  Path to the host FASTA genome file
  --host_bwa_index             [string]  Path to host genome directory or tar.gz archive for pre-built BWA index.
  --host_bowtie2_index         [string]  Path to host genome directory or tar.gz archive for pre-built BOWTIE2 index.
  --save_reference             [boolean] If generated by the pipeline save the BWA index in the results directory.

Read trimming options
  --trimmer                    [string]  Specifies the alignment algorithm to use - available options are 'fastp', 'trimmomatic'. [default: fastp]
  --adapters                   [string]  Path to FASTA adapters file
  --leading                    [integer] Instructs Trimmomatic to cut bases off the start of a read, if below a threshold quality [default: 3]
  --trailing                   [integer] Instructs Trimmomatic to cut bases off the end of a read, if below a threshold quality [default: 3]
  --average_quality            [integer] Instructs Trimmomatic or Fastp the average quality required in the sliding window [default: 20]
  --min_length                 [integer] Instructs Trimmomatic or Fastp to drop the read if it is below a specified length [default: 20]
  --qualified_quality_phred    [integer] Instructs Fastp to apply the --qualified_quality_phred option [default: 30]
  --unqualified_percent_limit  [integer] Instructs Fastp to apply the --unqualified_percent_limit option [default: 10]
  --skip_trimming              [boolean] Skip the adapter trimming step.
  --save_trimmed               [boolean] Save the trimmed FastQ files in the results directory.
  --save_trimmed_fail          [boolean] Save failed trimmed reads.

Alignment options
  --aligner                    [string]  Specifies the alignment algorithm to use - available options are 'bwa', 'bowtie2'. [default: bwa]
  --seq_center                 [string]  Sequencing center information to be added to read group of BAM files.
  --save_align_intermeds       [boolean] Save the intermediate BAM files from the alignment step.
  --skip_markduplicates        [boolean] Skip picard MarkDuplicates step.
  --skip_alignment             [boolean] Skip all of the alignment-based processes within the pipeline.
  --min_mapped                 [integer] Minimum number of mapped reads to be used as threshold to drop low mapped samples [default: 200]

Amplicon trimming options
  --primer_scheme_version      [string]  PrimalScheme RVFV primer scheme to use 'V1', 'V2' and 'V3'
  --ivar_trim_noprimer         [boolean] Unset -e parameter for ivar trim. Reads with primers are excluded by default
  --ivar_trim_min_len          [integer] Minimum length of read to retain after trimming [default: 20]
  --ivar_trim_min_qual         [integer] Minimum quality threshold for sliding window to pass [default: 20]
  --ivar_trim_window_width     [integer] Size of the sliding window [default: 4]
  --amplicon_left_suffix       [string]  Left suffix string in the amplicons primer bed file [default: _LEFT]
  --amplicon_right_suffix      [string]  Right suffix string in the amplicons primer bed file [default: _RIGHT]

Variant calling options
  --mpileup_depth              [integer] SAMtools mpileup max per-file depth, avoids excessive memory usage
  --min_base_quality           [integer] Skip bases with baseQ/BAQ smaller than this value when performing variant calling [default: 20]
  --min_coverage               [integer] Skip positions with an overall read depth smaller than this value when performing variant calling [default: 10]
  --min_allele_freq            [number]  Minimum allele frequency threshold for calling variant  [default: 0.25]
  --max_allele_freq            [number]  Maximum allele frequency threshold for calling variant  [default: 0.75]
  --save_mplieup               [boolean] Save SAMtools mpileup output file

Process skipping options
  --skip_multiqc               [boolean] Skip MultiQC.
  --skip_qc                    [boolean] Skip all QC steps except for MultiQC.
```

## Output

All the output results will be written to the `results` directory if no `--outdir` is not used. Masked and non-masked consensus genomes will be located in `bcftools/consensus`


See [usage docs](https://github.com/ajodeh-juma/rvfvampliconseq/usage) for all of the available options when running the pipeline.

## Pipeline Summary

By default, the pipeline currently performs the following:

* Sequencing quality control (`FastQC`)
* Quality control and preprocessing (`fastp`) or (`trimmomatic`)
* Reads alignment/mapping (`BWA`) or (`Bowtie2`)
* Alignment summary (`SAMtools`)
* Call variants (`iVar`)
* Annotate variants (`SnpEff`) or (`SnpSift`)
* Genome coverage (`BEDTools`)
* Visualization (`R`), (`ggplot2`)
* Overall pipeline run summaries (`MultiQC`)

## Documentation
Generating whole genome sequences of segmented viruses has largely depended on sequencing of partial gene sequences of the viruses. Here we implement a pipeline
that can be adopted to other segmented viruses in order to assemble complete genomic sequences from RNA metagenomic sequencing. We implement this pipeline to
generate complete genome sequences of Rift Valley fever virus, a tripartite virus having 3 segments - Small (S), Medium (M) and Large (L).

The pipeline comes bundled with reference genome and annotation, and the user only has to specify the segment to obtain 
full genome sequences. The pipeline calls variants
using `iVar` and annotates the variants using `SnpEff` and `SnpSift`

## Credits

**rvfvampliconseq** was originally written by @ajodeh-juma with inspiration from the @nf-core team, particularly on viralrecon

We thank the following people for their extensive assistance in the development
of this pipeline:

## License
rvfvampliconseq is free software, licensed under [MIT](https://github.com/ajodeh-juma/rvfvampliconseq/blob/master/LICENSE).

## Issues
Please report any issues to the [issues page](https://github.com/ajodeh-juma/rvfvampliconseq/issues).

## Contribute
If you wish to fix a bug or add new features to the software we welcome Pull Requests. We use
[GitHub Flow style development](https://guides.github.com/introduction/flow/). Please fork the repo, make the change, then submit a Pull Request against out master branch, with details about what the change is and what it fixes/adds. 
We will then review your changes and merge them, or provide feedback on enhancements.


## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi. -->
<!-- If you use  nf-core/rvfvampliconseq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

In addition, references of tools and data used in this pipeline are as follows:
> Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data. http://www.bioinformatics.babraham.ac.uk/projects/fastqc
>
> Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. 
> _Bioinformatics_, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560
>
> Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. 
> _Bioinformatics (Oxford, England)_, 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170
>
> Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. 
> _Bioinformatics_, 25(14), 1754–1760. https://doi.org/10.1093/bioinformatics/btp324
>
> Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. _Nat Methods_. 2012;9(4):357-359. Published 2012 Mar 4. doi:10.1038/nmeth.1923
>
> Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., & 1000 Genome Project Data Processing Subgroup. (2009). 
> The Sequence Alignment/Map format and SAMtools. 
> _Bioinformatics_, 25(16)> 2078–2079. https://doi.org/10.1093/bioinformatics/btp352
>
> Li H. A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. 
> _Bioinformatics_. 2011;27(21):2987-2993. doi:10.1093/bioinformatics/btr509
>
> R Core Team. (2017). R: A language and environment for statistical computing. _R Foundation for Statistical Computing_,. https://www.R-project.org/
>
> Grubaugh, N.D.; Gangavarapu, K.; Quick, J.; Matteson, N.L.; De Jesus, J.G.; Main, B.J.; Tan, A.L.; Paul, L.M.; Brackney, D.E.; Grewal, S.; et al. An Amplicon-Based Sequencing Framework for Accurately Measuring Intrahost Virus Diversity Using PrimalSeq and IVar. Genome Biol. 2019, 20, 8, doi:10.1186/s13059-018-1618-7
>
> Cingolani, P., Platts, A., Wang, l., Coon, M., Nguyen, T., Wang, L., Land, S. J., Lu, X., & Ruden, D. M. (2012). 
> A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. 
> _Fly_, 6(2), 80–92. https://doi.org/10.4161/fly.19695