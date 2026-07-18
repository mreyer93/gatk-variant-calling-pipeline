## Sequence read processing and alignment
*If you are starting with aligned bam files, skip this step!*
The sequence read processing pipeline performs the following steps: 
 - Fastqc and multiqc on raw data
 - Read quality trimming with `trim_galore`
 - Alignment to the reference genome with `bwa` and removal of duplicated reads
 - Various quality, coverage and alignment checks on the resulting bam files

First, copy the config file `process_reads/config_processing.yaml` to your working directory and define the run parameters. The `sample_file` defines sample names and paired-end fastq files. 

Next, gather reference files: the GRCh38 human genome fasta and a BED file of your project's capture-panel targets (there is no default target file - see [manual/requirements.md](requirements.md)).

You can launch the pipeline with a command like so:
```
snakemake -s /PATH/TO/pipeline/process_reads/processing.snakefile --configfile config_processing.yaml --use-conda --jobs 32 --cores 32 -k 
```
The `--use-conda` directive specifies that specific versions of the software used in this pipeline will be downloaded for each rule that requires them. To view a summary of which steps will be executed before running the pipeline, add the `-n` flag for a dry-run.

### Output files
The following folders and files are produced: 
```
00_qc_reports:
  pre_fastqc: fastqc reports on raw reads
  pre_multiqc: multiqc report on all samples, raw reads
  post_fastqc: fastqc reports on trimmed reads
  post_multiqc: multiqc report on all samples, trimmed reads
00_read_symlinks: symlinks to raw data for easy access
01_trimmed: sequencing reads processed with trim_galore
02_align:
  bam: reads aligned to the reference genome, including a version with duplicates removed
  cov: coverage statistics of specified genome intervals
  offtarget: bamfiles for reads that aligned outside of the specified intervals
  aligned_counts.txt: count of aligned reads per sample
  flagstat_offtarget.txt: samtools flagstat run on the offtarget bam files
primer_check.pdf: coverage statistics from the alignments
```