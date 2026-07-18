'''
GATK germline variant calling pipeline

This pipeline bam files aligend to the human reference genome and
performs the following steps in the GATK best practices workflow: 
 - Adds read groups and changes chromosome annotations to 'chr1' format, to normalize all input files
 - Base Qualtiy Score Recalibration (BQSR)
 - Germline variant calling using GATK haplotypecaller
 - genomicsDBimport
 - joint genotyping
 - variant filtration
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

# list of chromosomes to process 
chromosome_list = [f"chr{i}" for i in list(range(1,23))]

# source scripts for various functions
# setup, read config and sample file
include: "scripts/setup_germline.smk"
# prepare the bamfiles, including read groups and chr names
include: "scripts/prepare_bam_germline.smk"
# variant calling steps: haplotypecaller
include: "scripts/variant_calling_germline.smk"
# variant annotation steps: GATK Funcotator
include: "scripts/variant_annotation_germline.smk"

# add output files desired if not annotation_only setting
variant_files = []
if not config['annotation_only']:
    if not config['final_bam']:
        if config['skip_bqsr']:
            variant_files.append(expand(join(outdir, '01_prepare_bam/add_readgroups/{sample}.fixchr.bam'), sample=sample_list))
        else:
            variant_files.append(expand(join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'), sample=sample_list))
    variant_files.append(expand(join(outdir, '05_haplotypecaller/{sample}.g.vcf.gz'), sample=sample_list))
    variant_files.append(expand(join(outdir, '06_GDB/{chromosome}/FILE'), chromosome=chromosome_list))
    variant_files.append(expand(join(outdir, '06_GDB/{chromosome}/{chromosome}.vcf'), chromosome=chromosome_list))
    variant_files.append(expand(join(outdir, '07_roh/{chromosome}/{sample}_roh.txt.gz'), chromosome=chromosome_list, sample=sample_list))
    variant_files.append(expand(join(outdir, '08_roh_stats/{chromosome}_roh_stats.tsv'), chromosome=chromosome_list))
    
annotation_files=[]
if not config['skip_annotation']:
    annotation_files.append(join(outdir, '07_joint_vcf/02_variant_annotations/annotations_combined.vcf'))
    # annotation_files.append(join(outdir, '08_germline_metrics/done.tmp'))

rule all:
    input:
        variant_files,
        annotation_files
