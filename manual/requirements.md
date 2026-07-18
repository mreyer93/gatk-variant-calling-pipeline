### Requirements
The only software you have to install directly is Snakemake; everything else the pipeline calls
(GATK, samtools, bwa, R packages, etc.) is installed into a conda environment. It's easiest to do
this with [Conda](https://docs.conda.io/en/latest/miniconda.html) - install it by following the
instructions from that link.

I also recommend installing `mamba` for an improved user experience with conda. After conda is
installed and configured, run this command in your base environment.
```
conda install mamba -n base -c conda-forge
```

Then, create an environment named `gatk-pipeline` to use the variant calling pipeline. This will
use the specifications defined in the `envs/environment.yml` file.
```
mamba env create -f envs/environment.yml -n gatk-pipeline
```

You have to activate this environment whenever you want to run the pipeline.
```
conda activate gatk-pipeline
```

You also need the reference databases used by this pipeline available on your machine. The
pipeline targets **GRCh38**. None of these are bundled with this repo - download them per
project:
 - Human genome reference fasta, GRCh38, contigs named `chr1` etc. (the standard analysis-set
   fasta already uses this naming). Broad resource bundle:
   `gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta`
 - A BED file of your project's capture-panel targets (or a whole-genome/whole-exome target file -
   see the comment in `config_call_bam_GATK.yaml` for how to generate one) - there is no default,
   since a different project may use a different panel.
 - dbSNP and gnomAD resources - see the comments in `call_bam_GATK/config_call_bam_GATK.yaml` for
   the specific files and where to get them (Broad's public GATK resource bundles).
 - Funcotator data sources: run `gatk FuncotatorDataSourceDownloader --germline --hg38
   --validate-integrity --extract-after-download` (add `--somatic` too if you want the somatic
   data sources) rather than hardcoding a path - GATK 4.6.2 moved the hosting location for these,
   so a pinned URL/path is likely to go stale.

#### Install CADD scripts and download databases
CADD (v1.6+, which supports GRCh38) is used for variant deleteriousness scoring. The databases are
very large (400G+) and take time to download. Once the `gatk-pipeline` conda environment above is
created and activated, you will have access to the `cadd.sh` and `cadd-install.sh` scripts. Run
`cadd-install.sh` and install the databases you need (answer "no" to the first question, which
asks about installing a separate CADD conda environment - you already have one). You will need to
do this before running the pipeline all the way through.
