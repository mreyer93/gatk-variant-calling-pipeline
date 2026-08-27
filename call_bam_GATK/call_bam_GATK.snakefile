'''
GATK somatic variant calling pipeline

This pipeline bam files aligend to the human reference genome and
performs the following steps in the GATK best practices workflow: 
 - Adds read groups and changes chromosome annotations to 'chr1' format, to normalize all input files
 - Base Qualtiy Score Recalibration (BQSR)
 - Sample vs reference calling with Mutect2 for all samples
 - Tumor vs Normal calling with Mutect2 for all patient-timepoints that have paired samples
 - variant annotation for each of the "tumor vs reference" and "tumor vs normal" sets, including
    - CADD scoring
    - GATK VariantAnnotator 
    - GATK Funcotator
 - Variant annotations are combined into a single file with the most relevant columns kept. 
 - Variants in the "tumor vs reference" set are scored to bring the most relevant variants to the top of the list. 
 - The FACETS algorithm is used to call large copy number alterations and loss of heterozygosity events
'''

###############################################################################
import fnmatch as fn
import pandas as pd
from os.path import join, splitext
from collections import defaultdict
from snakemake.io import Wildcards, expand
import snakemake.rules
import sys
################################################################################

# source scripts for various functions
# setup, read config and sample file
include: "scripts/setup.smk"
# prepare the bamfiles, including read groups and chr names
include: "scripts/prepare_bam.smk"
# variant calling steps: tumor vs reference
include: "scripts/variant_calling_TvR.smk"
# variant calling steps: tumor vs normal
include: "scripts/variant_calling_TvN.smk"
# annotation steps: tumor vs reference
include: "scripts/variant_annotation_TvR.smk"
# annotation steps: tumor vs normal
include: "scripts/variant_annotation_TvN.smk"
# FACETS plotting CNA and LOH
include: "scripts/FACETS.smk"
# client-facing summary report (HTML/PDF) over whatever this run produced
include: "scripts/report.smk"


# add output files desired if not annotation_only setting
variant_files = []
if not config['annotation_only']:
    if config['skip_bqsr']:
        variant_files.append(expand(join(outdir, '01_prepare_bam/add_readgroups/{sample}.fixchr.bam'), sample=sample_list))
    else:
        variant_files.append(expand(join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'), sample=sample_list))
    variant_files.append(expand(join(outdir, '01_prepare_bam/depth_smoothed/{sample}.tsv'), sample=sample_list))
    variant_files.append(join(outdir, '01_prepare_bam/depth_smoothed_join.tsv'))
    variant_files.append(expand(join(outdir, '01_prepare_bam/GetPileupSummaries/{sample}_GetPileupSummaries.table'), sample=sample_list))
    variant_files.append(expand(join(outdir, '01_prepare_bam/CalculateContamination/{t_sample}_contam.table'), t_sample=all_t_contamination_samples))
    variant_files.append(expand(join(outdir, '02_variants_reference/01_gatk_variant_calling/filtered/{sample}_filtered.vcf'), sample=sample_list))
    #variant_files.append(expand(join(outdir, '03_variants_TvN/01_gatk_variant_calling/filtered/{patient_tp}_filtered.vcf'), patient_tp=tn_pt_tps))

annotation_files = []
if not config['skip_annotation']:
        annotation_files.append(expand(join(outdir, '02_variants_reference/02_variant_annotations/annotate_CADD/simple_table/{sample}_filtered_CADD_simple.vcf'), sample=sample_list))
        annotation_files.append(expand(join(outdir, '02_variants_reference/02_variant_annotations/annotate_Funcotator/{sample}_filtered_func.vcf'), sample=sample_list))
        annotation_files.append(expand(join(outdir, '02_variants_reference/02_variant_annotations/annotations_combined/filtered/{sample}_filtered_annotated.vcf'), sample=sample_list))
        annotation_files.append(expand(join(outdir, '02_variants_reference/03_variant_new_scores/prescore/filtered/{pt}.vcf'), pt=unique_pts))
        annotation_files.append(expand(join(outdir, '02_variants_reference/03_variant_new_scores/prescore/unfiltered/{pt}.vcf'), pt=unique_pts))
        annotation_files.append(expand(join(outdir, '02_variants_reference/03_variant_new_scores/final/filtered/{pt}.vcf'), pt=unique_pts))
        # annotation_files.append(expand(join(outdir, '03_variants_TvN/02_variant_annotations/annotations_combined/filtered/{patient_tp}_filtered_annotated.vcf'), patient_tp=tn_pt_tps))

report_files = []
if config['make_report']:
    report_files.append(join(outdir, '09_report/somatic_report.html'))
    if config['report_pdf']:
        report_files.append(join(outdir, '09_report/somatic_report.pdf'))

rule all:
    input:
        variant_files,
        annotation_files,
        report_files,
        #expand(join(outdir, '04_FACETS/{patient_tp}.pdf'), patient_tp=tn_pt_tps)
