# Variant calling with GATK

A reproducible, scalable Snakemake pipeline for calling somatic and germline variants in
targeted/exome sequencing data using the Genome Analysis Toolkit (GATK), plus a companion
pipeline for QC processing of raw sequencing reads and alignment to a reference genome.

This is a general-purpose pipeline meant to be reused across projects: bring your own reference
genome build, capture panel, and sample sheet via the config file rather than relying on any
built-in defaults.

## Quickstart
Install miniconda and mamba. Create a conda environment from the file in `envs/environment.yml`:
```
mamba env create -f envs/environment.yml -n gatk-pipeline
conda activate gatk-pipeline
```
Download the required references (see [manual/requirements.md](manual/requirements.md) - the
pipeline targets GRCh38). Copy the config file at `call_bam_GATK/config_call_bam_GATK.yaml` to
your working directory and change the settings to fit your desired output directory, metadata
files, capture-panel targets, and reference locations. Then launch the pipeline, changing the
resource allocations to match what you have available on your machine:
```
snakemake -s /PATH/TO/pipeline/call_bam_GATK/call_bam_GATK.snakefile --configfile config_call_bam_GATK.yaml --jobs 32 --cores 32 -k
```

A similar process runs germline calling (`call_bam_GATK_germline.snakefile` +
`config_call_bam_GATK_germline.yaml`) and the sequence read QC/alignment pipeline, by copying the
config file at `process_reads/config_processing.yaml` to your working directory, changing the
parameters, and running:
```
snakemake -s /PATH/TO/pipeline/process_reads/processing.snakefile --configfile config_processing.yaml --jobs 32 --cores 32 -k
```

Don't have 400G+ free for CADD? Use `config_call_bam_GATK_local.yaml` /
`config_call_bam_GATK_germline_local.yaml` instead - see [manual/requirements.md](manual/requirements.md)
for what that trades off, or [cloud/gcp/README.md](cloud/gcp/README.md) to run the full
pipeline on a GCP VM instead of locally.

## User Manual
1. [Installation requirements](manual/requirements.md)
2. [Sequence read QC and alignment](manual/processing.md)
3. [Somatic variant calling with GATK](manual/somatic.md)
4. [Germline variant calling with GATK](manual/germline.md)
5. [Other things in this repository](manual/other.md)

## Repository layout
- `call_bam_GATK/` - somatic and germline variant-calling pipelines (the main pipelines here)
- `process_reads/` - read QC, trimming, and alignment pipeline
- `envs/` - conda/mamba environment definitions
- `manual/` - detailed docs for each pipeline stage
- `cloud/gcp/` - scripts to run the full pipeline on a GCP VM (see [cloud/gcp/README.md](cloud/gcp/README.md))
- `archive/` - historical, project-specific analyses and superseded code, kept for reference only
  (not maintained, not guaranteed to run) - see [archive/README.md](archive/README.md)
