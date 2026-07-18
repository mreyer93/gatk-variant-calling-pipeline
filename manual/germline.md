### Germline variant calling
Germline variant calling uses GATK HaplotypeCaller, joint-genotyped across all samples per
chromosome (via GenomicsDBImport), which is the current GATK best-practices approach and keeps
computation manageable on large batches. It operates on bam files (either from the processing
pipeline or another source) and performs the following steps:
1. Adds read groups and changes chromosome annotations to 'chr1' format, to normalize all input files
2. Base Quality Score Recalibration (BQSR)
3. Germline variant calling per-sample with HaplotypeCaller (GVCF mode)
4. GenomicsDBImport and joint genotyping, run per-chromosome
5. Hard-filtering of the joint-genotyped SNPs
6. Runs of homozygosity (ROH) detection with `bcftools roh`, per chromosome per sample, plus
   aggregated ROH stats
7. Variant annotation with CADD scoring and GATK Funcotator, combined into a single annotated file

#### Configure the pipeline
Copy `call_bam_GATK/config_call_bam_GATK_germline.yaml` to your working directory and change the parameters to fit your run. The `bam_metadata` file only needs two columns:

sample, bamfile

Reference files and the capture-panel target BED must be specified - see [manual/requirements.md](requirements.md) (there is no default target file, since different projects use different panels). The CADD scoring install from that page must also be performed. If you have already created the variant files and just want to re-run annotation with different parameters, the `annotation_only` flag in the config file allows for this.

Optional config flags (all default to `False` if omitted): `skip_bqsr` (use the read-group-fixed bam directly, skipping BQSR), `final_bam` (use the input bam as-is, skipping both BQSR and chromosome-fixing), `skip_annotation` (stop after variant calling, skip the annotation step).

#### Run the pipeline
```
snakemake -s /PATH/TO/pipeline/call_bam_GATK/call_bam_GATK_germline.snakefile --configfile config_call_bam_GATK_germline.yaml --use-conda --jobs 32 --cores 32 -k
```

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
  germline_calls_SNP_hard_filter_select.vcf: hard-filtered, SNP-only joint-genotyped calls
  02_variant_annotations
    annotations_combined.vcf: CADD + Funcotator annotations combined into one file
```
