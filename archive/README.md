# Archive

This directory holds everything from the original repo that is **not** part of the generalized,
reusable pipeline (see the top-level [README.md](../README.md) for that). Nothing in here is
maintained, documented, or guaranteed to run — it's kept only as historical reference and for
occasional code-reuse when tackling a similar analysis.

- `p6_investigation/` — deep-dive into a specific patient/cohort haplotype and ROH comparisons
  against 1000 Genomes. Hardcoded to specific sample IDs.
- `EGA_recalling/` — recalling variants in a specific EGA cohort (Sentieon + gnomAD queries via
  Hail).
- `EGA_variant_processing/` — germline/somatic variant annotation and LOH association analysis
  for that same EGA cohort.
- `rescreen_gnomad_comparison/` — one notebook comparing a specific rescreen dataset to gnomAD.
- `rnaseq/` — DESeq2/MSigDB/KEGG analysis for two named mouse RNA-seq experiments.
- `splicing/` — GFI1/SF3B1 splicing investigation in the EGA data.
- `vaf_decrease_plotting/` — figure generation for VAF/LOH-decrease-over-time plots for a specific
  report.
- `parse_bio/` — single ad hoc script, unclear scope.
- `scripts_misc/` — sample-specific one-off scripts that used to live in the top-level `scripts/`
  directory (coverage checks for particular samples/genes).
- `exploratory_analysis/` — ad hoc, hardcoded-sample-ID scratch scripts for LOH calling and VCF
  clustering. These represent ideas for pipeline features that were never generalized (see
  `manual/other.md` in the main pipeline for the original design note on VCF clustering) —
  worth revisiting as a properly parameterized Snakemake step if needed again, rather than reusing
  as-is.
- `deprecated_aws/` — the old S3/EC2-based execution model for `call_bam_GATK`, superseded by the
  current local-BAM Snakemake pipeline. A leaked AWS access key and an EGA login that were
  originally in `aws_config.sh` have been removed (rotate/deactivate that AWS key if it's still
  live). AWS/cloud execution support will be redesigned from scratch as a separate effort rather
  than resurrecting this.
