### Somatic variant calling 
The somatic variant calling is the best developed section of this pipeline, and operates on bam files (either from the processing pipeline or another source). The pipeline implements the [GATK best practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-) and performs the following steps: 
1. Adds read groups and changes chromosome annotations to 'chr1' format, to normalize all input files
2. Base Qualtiy Score Recalibration (BQSR)
3. Sample vs reference calling with Mutect2 for all samples
4. Tumor vs Normal calling with Mutect2 for all patient-timepoints that have paired samples
5. variant annotation for each of the "tumor vs reference" and "tumor vs normal" sets, including
    - CADD scoring
    - GATK VariantAnnotator 
    - GATK Funcotator
6. Variant annotations are combined into a single file with the most relevant columns kept. 
7. Variants in the "tumor vs reference" set are scored to bring the most relevant variants to the top of the list. 
8. Calling large copy number alterations and loss of heterozygosity events with the [FACETS](https://github.com/mskcc/facets/) algorithm.

#### Configure the pipeline
Copy `call_bam_GATK/config_call_bam_GATK.yaml` to your working directory and change the parameters to fit your run. The `bam_metadata` file gives all the input files and annotations. This must be a tab-delimited file with the following column headings:

sample, bamfile, patient, timepoint, tumor_normal, library, platform, platform_unit

The last three columns are required for adding read groups to the bam files, but not used in any other calculations. I set them all to "lib1    ILLUMINA    unit1"

Reference files, the capture-panel target BED, and annotation directories must be specified - see [manual/requirements.md](requirements.md) for where to get each one (there is no default target file, since different projects use different panels). The CADD scoring install from that page must also be performed. If you have already created the variant files and just want to run the annotation again, perhaps with different parameters, the `annotation_only` flag in the config file allows for this.

Mutect2 is run with `--germline-resource` pointed at the gnomAD resource in the config, to help distinguish somatic from germline calls, and orientation-bias correction (`LearnReadOrientationModel`/`--ob-priors`) is applied on both the tumor-vs-reference and tumor-vs-normal paths. There is no panel-of-normals (PoN) support yet - if you have a cohort of matched normal samples, building and wiring in a PoN is a worthwhile follow-on improvement for more accurate somatic filtering.

##### Running without CADD (e.g. on a laptop)
CADD's database is 400G+, which won't fit on most local machines. Unlike the germline
pipeline, there's no lightweight substitute here: the "new score" ranking below multiplies
CADD's score by AF and fold-change, and the `filtered` output requires a CADD score to be
present at all. AlphaMissense (used as CADD's local substitute in the germline pipeline) only
scores missense SNVs, so swapping it in here would silently drop every indel, nonsense, and
splice-site call from the scored/filtered output - often the clinically important calls in a
cancer gene panel. That's worse than just not scoring locally. Set `skip_annotation: True`
instead (`config_call_bam_GATK_local.yaml` is set up this way) to get variant calls, filtering,
and FACETS without CADD/Funcotator/scoring, and run the full pipeline with CADD on a machine
that has the disk space for it.

#### Run the pipeline
Similar to the other pipelines, you can run somatic variant calling with a command like:
```
snakemake -s /PATH/TO/pipeline/call_bam_GATK/call_bam_GATK.snakefile --configfile config_call_bam_GATK.yaml --use-conda --jobs 32 --cores 32 -k 
```
Lower `--jobs`/`--cores` to match your machine - e.g. `--jobs 4 --cores 4` on a 4-core laptop.
The `--use-conda` directive specifies that specific versions of the software used in this pipeline will be downloaded for each rule that requires them. To view a summary of which steps will be executed before running the pipeline, add the `-n` flag for a dry-run.

#### New score files
To pull out variants that are likely to be true meaningful somatic events, I developed the following methodology. This is run for each patient-timepoint that has at least one T and one N sample.
- In cases where a variant appears in a T sample but not a N sample, force a genotype call at that locus in the normal sample. This will differentiate if the variant was not called in the N sample because of low depth, or it was truly not present. 
- Rank variants by the following metric (called the _new score_)
    - Calculate the log2 fold change of allelic fraction from N to T (plus epsilon to account for zero AF in N samples). The mean of samples in each condition is used. 
    - _new score_  = CADD score * mean Tumor AF * log2 fold change

Files in the `03_variant_new_scores` directory are processed with this methodology and should pull high-confidence variants to the top of the list. 

#### Client-facing summary report
The pipeline finishes by building a self-contained report summarising the whole run, intended to be handed directly to a collaborator or client:
```
09_report/somatic_report.html    # self-contained, opens in any browser
09_report/somatic_report.pdf     # same content, paginated for formal delivery
09_report/summary_tables/*.tsv   # every figure's underlying table, for Excel/re-analysis
```
It covers: run configuration, sequencing depth and tumour-in-normal contamination QC, variant counts per sample, functional-consequence and variant-type breakdowns, allele-fraction and CADD distributions, recurrently mutated genes, a ranked shortlist of prioritised variants, and FACETS purity — plus a methods/provenance section and explicit interpretation caveats.

The report adapts to how the run was configured. With `skip_annotation: True` (the local config), the gene-level, consequence, and deleteriousness sections are omitted automatically and the QC and variant-count sections still render — it does not fail. Sections whose inputs are missing for any other reason (no matched normals, FACETS not run) are likewise skipped with a short note rather than erroring.

Control it with these config options (all optional):
```yaml
make_report: True                       # build the report at all
report_pdf: True                        # also render PDF (needs LaTeX/tectonic from envs/rmarkdown.yml)
project_name: "Client X - AML panel"    # report title
report_top_n: 40                        # rows in the prioritised-variant table
```
HTML and PDF are rendered by separate Snakemake rules, so a broken LaTeX toolchain costs you the PDF but never the HTML. Set `report_pdf: False` if you don't need PDF.

#### Output files
The following output directories and files are generated from this pipeline: 
```
01_prepare_bam
  add_readgroups: bam files with readgroups added.
  CalculateContamination: results from GATK CalculateContamination for calculating tumor/normal contamination
  GetPileupSummaries: results from GATK GetPileupSummaries
  recalibrate: results of Base Quality Score Recalibration (BQSR), used in the next steps

02_variants_reference: results of variant calling vs the fasta reference file
  01_gatk_variant_calling:
    filtered: variant calling results with GATK filters calculated. Variants failing the filters remain in the file, with an annotation in the FILTER field. 
    unfiltered: raw variant calling results 
  
  02_variant_annotations
    annotate_CADD: CADD scoring run on each variant file
    annotate_Funcotator: GATK Funcotator run on each variant file
    annotations_combined: 
      filtered: A combined annotated vcf-like file for each sample. Filters applied to remove low-quality variants
      unfiltered: same, but without the filters applied

  03_variant_new_scores: variant results processed with the new score methodolgy (see above)
    final: 
      filtered: new scored variants, filtered to min AF and depth as specified in the config file
      filtered_GATK: Same as above, but variants must pass the GATK filters listed in the config file
    prescore: processing files, including initial scoreing and genotyping of normals with missing variant calls

03_variants_TvN: Result of running GATK Mutect2 on paired tumor and normal bam files. One output is produced for each patient-timepoint. These calls are generally less reliable unless the Normal samples are very clean (due to tumor-in-normal contamination).
  01_gatk_variant_calling
    filtered: variant calling results with GATK filters calculated. Variants failing the filters remain in the file, with an annotation in the FILTER field. 
    unfiltered: raw variant calling results 
  02_variant_annotations
    annotate_CADD: CADD scoring run on each variant file
    annotate_Funcotator: GATK Funcotator run on each variant file
    annotations_combined: 
      filtered: A combined annotated vcf-like file for each sample. Filters applied to remove low-quality variants
      unfiltered: same, but without the filters applied

04_FACETS: data file and plots from running the FACETS algorithm for LOH and CNA detection.

09_report: client-facing summary report (see above)
  somatic_report.html / .pdf: the report itself
  summary_tables: TSV table behind every figure in the report
  logs: rendering logs, useful if a report build fails
```