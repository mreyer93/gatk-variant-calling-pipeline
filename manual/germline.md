### Germline variant calling
Germline variant calling uses GATK HaplotypeCaller, joint-genotyped across all samples per
chromosome (via GenomicsDBImport), which is the current GATK best-practices approach and keeps
computation manageable on large batches. It operates on bam files (either from the processing
pipeline or another source) and performs the following steps:
1. Adds read groups and changes chromosome annotations to 'chr1' format, to normalize all input files
2. Base Quality Score Recalibration (BQSR)
3. Germline variant calling per-sample with HaplotypeCaller (GVCF mode)
4. GenomicsDBImport and joint genotyping, run per-chromosome
5. Hard-filtering of the joint-genotyped SNPs and indels, each with GATK's recommended
   thresholds ([Hard-filtering germline short variants](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112))
6. Runs of homozygosity (ROH) detection with `bcftools roh`, per chromosome per sample, plus
   aggregated ROH stats
7. Variant annotation with CADD scoring (or AlphaMissense - see below) and GATK Funcotator,
   combined into a single annotated file

#### Configure the pipeline
Copy `call_bam_GATK/config_call_bam_GATK_germline.yaml` to your working directory and change the parameters to fit your run. The `bam_metadata` file only needs two columns:

sample, bamfile

Reference files and the capture-panel target BED must be specified - see [manual/requirements.md](requirements.md) (there is no default target file, since different projects use different panels). The CADD scoring install from that page must also be performed. If you have already created the variant files and just want to re-run annotation with different parameters, the `annotation_only` flag in the config file allows for this.

Optional config flags (all default to `False` if omitted): `skip_bqsr` (use the read-group-fixed bam directly, skipping BQSR), `final_bam` (use the input bam as-is, skipping both BQSR and chromosome-fixing; it must be indexed, as `.bam.bai` or `.bai`), `skip_annotation` (stop after the filtered calls, skip the annotation step).

`chromosomes` lists the chromosomes to joint-genotype, one GenomicsDBImport workspace each; it defaults to `chr1`-`chr22`. Every one listed must be in the reference. Adding `chrX`/`chrY` also needs ploidy handling for male samples (HaplotypeCaller's `-ploidy`), which the pipeline does not do yet.

##### CADD licensing: commercial and client work
CADD is free for non-commercial use only, and the full config runs it. For paid client work
without a CADD licence use `config_call_bam_GATK_germline_commercial.yaml`, which sets
`skip_cadd: True` and points at AlphaMissense. Unlike the somatic commercial config, this one
keeps the rest of the annotation stage: Funcotator still runs, and only the CADD step is
skipped. See [manual/requirements.md](requirements.md#annotation-licensing-and-commercial-use)
for the licence terms on both CADD and AlphaMissense.

##### Running without CADD (e.g. on a laptop)
CADD's database is 400G+, which won't fit on most local machines. Set `skip_cadd: True` and
point `alphamissense_file` at a tabix-indexed `AlphaMissense_hg38.tsv.gz`
(`./scripts/download_references.sh --with-alphamissense`, ~640M) to get
[AlphaMissense](https://www.nature.com/articles/s41586-023-06213-3) deleteriousness scores
instead - it's also the current [ClinGen SVI](https://clinicalgenome.org/)-preferred predictor
for germline missense variant classification (ACMG PP3/BP4), so this isn't just a fallback.
The catch: AlphaMissense only scores missense SNVs - indels, nonsense, splice-site, and
non-coding variants will have no `AM_pathogenicity`/`AM_class` value even with this enabled.
`config_call_bam_GATK_germline_local.yaml` is a ready-to-use config set up this way. This
tradeoff is specific to germline - see `manual/somatic.md` for why the somatic pipeline's
local config skips deleteriousness scoring entirely instead of doing the same substitution.

#### Test it
`./test/run_test_germline.sh` runs the whole pipeline (annotation off) on the small public
tumour/normal pair the somatic smoke test uses, treated as two germline samples, and checks
that the outputs are well formed. It does not measure accuracy: that data has no truth set.

#### Run the pipeline
```
snakemake -s /PATH/TO/pipeline/call_bam_GATK/call_bam_GATK_germline.snakefile --configfile config_call_bam_GATK_germline.yaml --use-conda --jobs 32 --cores 32 -k
```
Lower `--jobs`/`--cores` to match your machine - e.g. `--jobs 4 --cores 4` on a 4-core laptop.

#### Output files
```
01_prepare_bam
  add_readgroups: bam files with readgroups added and chromosome names normalized
  recalibrate: results of Base Quality Score Recalibration (BQSR), used in the next steps

05_haplotypecaller: per-sample GVCFs from HaplotypeCaller

06_GDB: per-chromosome GenomicsDBImport workspaces and joint-genotyped VCFs

07_roh: per-chromosome, per-sample runs-of-homozygosity calls (bcftools roh)
08_roh_stats: aggregated ROH statistics per chromosome

07_joint_vcf
  germline_calls_hard_filter_select.vcf: joint-genotyped SNPs and indels that pass the hard filters
  germline_calls_hard_filter.vcf: the same calls before selection, with the failed filter named in FILTER
  02_variant_annotations
    annotations_combined.vcf: CADD (or AlphaMissense, if skip_cadd is set) + Funcotator
      annotations combined into one file - look for the CADD_phred column, or
      AM_pathogenicity/AM_class if running the AlphaMissense-based local config
```
