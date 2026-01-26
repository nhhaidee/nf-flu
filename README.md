# CFIA-NCFAD/nf-flu - Influenza A and B Virus Genome Assembly Nextflow Workflow

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15093332.svg)](https://doi.org/10.5281/zenodo.15093332)
[![CI](https://github.com/CFIA-NCFAD/nf-flu/actions/workflows/ci.yml/badge.svg)](https://github.com/CFIA-NCFAD/nf-flu/actions/workflows/ci.yml)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with podman](https://img.shields.io/badge/run%20with-podman-1d355c.svg?labelColor=000000)](https://podman.io/)

## Introduction

**nf-flu** is a [Nextflow][] bioinformatics analysis pipeline for assembly and analysis of Influenza A and B viruses from Illumina or Nanopore sequencing data or previously assembled FASTA sequences.
Since Influenza has a segmented genome consisting of 8 gene segments, the pipeline will automatically select the top matching reference sequence from NCBI for each gene segment based on [IRMA][] assembly and nucleotide [BLAST][] against all Influenza sequences from NCBI.
Users can also provide their own reference sequences to include in the top reference sequence selection process.
After reference sequence selection, the pipeline performs read mapping to each reference sequence, variant calling and depth-masked consensus sequence generation.

> **Note:** The officially supported version of the pipeline is [CFIA-NCFAD/nf-flu](https://github.com/CFIA-NCFAD/nf-flu). If you have issues with using the pipeline, please [create an issue](https://github.com/CFIA-NCFAD/nf-flu/issues/new/choose) on [CFIA-NCFAD/nf-flu](https://github.com/CFIA-NCFAD/nf-flu) repo.

## Pipeline summary

1. Download latest NCBI Influenza virus sequences and metadata (see [docs](docs/update_seqs_db.md) for more details).
2. Merge reads of re-sequenced samples ([`cat`](http://www.linfo.org/cat.html)) (if needed).
3. Assembly of Influenza gene segments with [IRMA][] using the built-in FLU module
4. Nucleotide [BLAST][] search against [NCBI Influenza DB][] sequences
5. H/N subtype prediction and Excel XLSX report generation based on BLAST results.
6. Automatically select top match reference sequences for segments
7. Read mapping, variant calling and consensus sequence generation for each segment against top reference sequence based on BLAST results.
8. Annotation of consensus sequences with [VADR][]
9. [FluMut][] detection of molecular markers and mutation in Influenza A(H5N1) viruses.
10. [GenoFLU][] genotyping of North American H5 viruses. [Genin2][] genotyping using information from clade 2.3.4.4b H5Nx viruses collected in Europe since October 2020.
11. HA cleavage site prediction and classification
12. [Nextclade][] clade assignment, mutation calling and sequence quality checks.
13. [MultiQC][] report generation.

![nf-flu workflow](assets/nf-flu-pipeline-diagram.svg)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`; latest stable release recommended!).
2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Apptainer`][], [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort)_
3. Download the pipeline and test it on a minimal dataset with a single command:

    For Illumina workflow test:

    ```bash
    nextflow run CFIA-NCFAD/nf-flu -profile test_illumina,<docker/apptainer/singularity/podman/shifter/charliecloud/conda> \
      --max_cpus $(nproc) # use all available CPUs; default is 2
    ```

    For Nanopore workflow test:

    ```bash
    nextflow run CFIA-NCFAD/nf-flu -profile test_nanopore,<docker/apptainer/singularity/podman/shifter/charliecloud/conda> \
      --max_cpus $(nproc) # use all available CPUs; default is 2
    ```

    > * If you are using `apptainer`/`singularity` then the pipeline will auto-detect this and attempt to download the Apptainer/Singularity images directly as opposed to performing a conversion from Docker images. If you are persistently observing issues downloading Apptainer/Singularity images directly due to timeout or network issues then please use the `--singularity_pull_docker_container` parameter to pull and convert the Docker image instead. Alternatively, it is highly recommended to use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to pre-download all of the required containers before running the pipeline and to set the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options to be able to store and re-use the images from a central location for future pipeline runs.
    > * If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Run your own analysis

    * [Optional] Generate an input samplesheet from a directory containing Illumina FASTQ files (e.g. `/path/to/illumina_run/Data/Intensities/Basecalls/`) with the included Python script [`fastq_dir_to_samplesheet.py`](bin/fastq_dir_to_samplesheet.py) **before** you run the pipeline (requires Python 3 installed locally) e.g.

        ```bash
        python ~/.nextflow/assets/CFIA-NCFAD/nf-flu/bin/fastq_dir_to_samplesheet.py \
          -i /path/to/illumina_run/Data/Intensities/Basecalls/ \
          -o samplesheet.csv
        ```

    * Typical command for Illumina sequencing data

        ```bash
        nextflow run CFIA-NCFAD/nf-flu \
          --input samplesheet.csv \
          --platform illumina \
          --profile <docker/apptainer/singularity/podman/shifter/charliecloud/conda>
        ```

    * Typical command for Nanopore sequencing data

      ```bash
      nextflow run CFIA-NCFAD/nf-flu \
        --input samplesheet.csv \
        --platform nanopore \
        --profile <docker/apptainer/singularity/conda>
      ```

    * Run analysis on FASTA files within a directory

      ```bash
      nextflow run CFIA-NCFAD/nf-flu \
        --input /path/to/fasta_files/ \
        --platform assemblies \
        --profile <docker/apptainer/singularity/conda>
      ```

## Documentation

The nf-flu pipeline comes with:

* [Usage](docs/usage.md) and
* [Output](docs/output.md) documentation.

## Note

* The pipeline requires internet access to download the latest NCBI Influenza virus sequences and metadata at the start of each run unless a local copy is provided via the `--ncbi_influenza_fasta` and `--ncbi_influenza_metadata` parameters. However, if you encounter error "Unable to stage foreign file" when online repository (Figshare) changes URL or applies security restriction for automated download, you can manually download DB (<https://api.figshare.com/v2/file/download/53449877> and <https://api.figshare.com/v2/file/download/53449874>) and store it locally and then provide the local paths to the downloaded files via the `--ncbi_influenza_fasta` and `--ncbi_influenza_metadata` parameters. If there is any problem, please open an issue.

## Contributors

* [Peter Kruczkiewicz](https://github.com/peterk87) ([CFIA-NCFAD](https://github.com/CFIA-NCFAD)) - lead developer
* [Hai Nguyen](https://github.com/nhhaidee) ([CFIA-NCFAD](https://github.com/CFIA-NCFAD)) - Nanopore workflow
* [Abdallah Meknas](https://github.com/ameknas-phac) (Influenza, Respiratory Viruses, and Coronavirus Section (IRVC), Public Health Agency of Canada (PHAC)) - expansion of the Illumina workflow
* [Cass Erdelyan](https://github.com/cerdelyan/) ([CFIA-NCFAD](https://github.com/CFIA-NCFAD)) - development, testing and valuable feedback

## Credits

* [nf-core](https://nf-co.re) project for establishing Nextflow workflow development best-practices, [nf-core tools](https://nf-co.re/tools-docs/) and [nf-core modules](https://github.com/nf-core/modules)
* [nf-core/viralrecon](https://github.com/nf-core/viralrecon) for inspiration and setting a high standard for viral sequence data analysis pipelines
* [Conda](https://docs.conda.io/projects/conda/en/latest/) and [Bioconda](https://bioconda.github.io/) project for making it easy to install, distribute and use bioinformatics software.
* [Biocontainers](https://biocontainers.pro/) for automatic creation of [Docker] and [Apptainer]/[Singularity] containers for bioinformatics software in [Bioconda]

[Apptainer]: https://apptainer.org/
[BcfTools]: https://samtools.github.io/bcftools/
[BLAST]: https://blast.ncbi.nlm.nih.gov/Blast.cgi
[Clair3]: https://github.com/HKU-BAL/Clair3
[Docker]: https://www.docker.com/
[FluMut]: https://github.com/izsvenezie-virology/FluMut
[Freebayes]: https://github.com/freebayes/freebayes
[Genin2]: https://github.com/izsvenezie-virology/genin2
[GenoFLU]: https://github.com/USDA-VS/GenoFLU
[IRMA]: https://wonder.cdc.gov/amd/flu/irma/
[Medaka]: https://github.com/nanoporetech/medaka
[Minimap2]: https://github.com/lh3/minimap2/
[Mosdepth]: https://github.com/brentp/mosdepth
[MultiQC]: https://multiqc.info/
[NCBI Influenza DB]: https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database
[NCBI Influenza Virus Resource]: https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database
[Nextclade]: https://clades.nextstrain.org/
[Nextflow]: https://www.nextflow.io/
[nf-core]: https://nf-co.re/
[Samtools]: https://www.htslib.org/
[seqtk]: https://github.com/lh3/seqtk
[Singularity]: https://www.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps
[table2asn]: https://www.ncbi.nlm.nih.gov/genbank/table2asn/
[VADR]: https://github.com/ncbi/vadr
