# worklfow implementing GATK best practices
# from an already aligned BAM file
# takes as input a mapping from sample to BAM file path
# and some information to be used in the future: 
    # patient, sample type, timepoint, seq lane, etc
###############################################################################
import fnmatch as fn
import pandas as pd
from os.path import join, splitext
from collections import defaultdict
from snakemake.io import Wildcards, expand
import snakemake.rules
import sys
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()
################################################################################
# define parameters from the configfile 
outdir = config['output_directory']
# for now, used the ref file with fixed chromosome names
REF_FILE = config['REF_FILE']
dbsnp_file =config['dbsnp_file']
# targets with 100bp windows for mutation calling
idt_target_w100 = config['idt_target_w100']
# genome version for references to use
# right now this is fixed with "GRCh37"
genome_version = "GRCh37"

metadata = pd.read_csv(config["bam_metadata"], delimiter = '\t')
req_columns = ['sample', 'EGAF', 'bamfile', 'patient', 'timepoint', 'tumor_normal', 'library', 'platform', 'platform_unit']
if not(all([r in list(metadata.columns) for r in req_columns])):
    sys.exit("Metadata must contain the following columns " + str(req_columns))
# set rownames to EGAF
metadata.index = list(metadata['EGAF'])

# ensure no duplicates in EGAF, bamfile
EGAF_list = list(metadata['EGAF'])
bamfile_list = list(metadata['bamfile'])
if (any(EGAF_list.count(x) > 1  for x in EGAF_list)):
    sys.exit("Check metadata, duplictes in EGAF column")
if (any(bamfile_list.count(x) > 1  for x in bamfile_list)):
    sys.exit("Check metadata, duplictes in bamfile column")
# final sample list is what we want to call variants on, this is all the unique samples
sample_list = list(set(metadata['sample']))
# list of each EGAF for each sample
sample_to_EGAF = {}
for s in sample_list:
    sample_to_EGAF[s] = list(metadata[metadata['sample']==s]['EGAF'])
print(sample_to_EGAF)


# print(sample_list)
# ensure all files specified actually exist
# might implement this with S3 in the future
# check_list = [not os.path.exists(a) for a in bamfile_list]
# if(any(check_list)):
#     missing_files = [b for a,b in zip(check_list, bamfile_list) if a]
#     print("Missing files: ")
#     for a in missing_files: print(a) 
#     sys.exit("Some specified bam files do not exist")

# map from EGAF ids to bamfiles
bam_map = {a: b for a,b in zip(EGAF_list, bamfile_list)}

# add annotation for patient-timepoint pair
metadata['pt_tp'] = [str(a)+'_'+str(b) for a, b in zip(metadata['patient'], metadata['timepoint'])]       
# sample name to final bamfile
sample_to_final_bam = {a:b for a, b in zip(list(metadata['sample']), expand(join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'), sample= list(metadata['sample'])))}
# map out what we need for calling tumor vs normal for each patient-timepoint pair
# this will be tumors from this timepoint, and the DX normals
tn_pt_tps = []
# map from patient-timepoint to samples and final bamfiles
pt_tp_to_samples = {}
pt_tp_to_final_bams = {}
# tumor samples for each patient-timepoint
t_samples_pt_tp = {}
# normal samples to use for each patient-timepoint, 
# this will be the same within a patient becasue normals 
# only collected at Dx
n_samples_pt_tp = {}
# also a version for all samples from a given patient
n_samples_pt = {}
t_samples_pt = {}


unique_pts = list(set(metadata['patient']))
for p in unique_pts:
    m_pt = metadata[metadata['patient']==p]
    pt_tps = list(set(m_pt['pt_tp']))
    if ('T' in list(m_pt['tumor_normal']) and 'N' in list(m_pt['tumor_normal'])):
        n_samples_pt[p] = list(m_pt[m_pt['tumor_normal']=='N']['sample'])
        t_samples_pt[p] = list(m_pt[m_pt['tumor_normal']=='T']['sample'])
        # iterate over all timepoints for this patient 
        for pt_tp in pt_tps:
            m_pt_tp = m_pt[m_pt['pt_tp']==pt_tp]
            tn_pt_tps += [pt_tp]
            # add all normal samples to this patient 
            n_samples_pt_tp[pt_tp] = list(m_pt[m_pt['tumor_normal']=='N']['sample'])
            t_samples_pt_tp[pt_tp] = list(m_pt_tp[m_pt_tp['tumor_normal']=='T']['sample'])
            # all samples to use in this calculation            
            pt_tp_to_samples[pt_tp] = n_samples_pt_tp[pt_tp] + t_samples_pt_tp[pt_tp]
            # all bams to use in this calculation
            pt_tp_to_final_bams[pt_tp] = expand(join('data-bucket/gatk_calling/', outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'), sample= pt_tp_to_samples[pt_tp])


# for tumor-normal contamination estimation, which pairs do we want to calculate
# it's each tumor vs one normal, ideally the best one
# for now just use the first? 
pt_contamination_normals = {a:b[0] for a,b in n_samples_pt.items()}
# this gets ruin for every T sample from each patient 
pt_contamination_tn_sets = {}
for pt,samples in t_samples_pt.items():
    for s in samples:
        pt_contamination_tn_sets[s] = [s, pt_contamination_normals[pt]]

all_t_contamination_samples = [item for sublist in list(t_samples_pt.values()) for item in sublist] 


rule all: 
    input: 
        # S3.remote(expand(join('data-bucket/gatk_calling', outdir, '02_variants_reference/unfiltered/{sample}_unfiltered.vcf'), sample=sample_list)),
        # S3.remote(expand(join('data-bucket/gatk_calling', outdir, '02_variants_reference/filtered/{sample}_filtered.vcf'), sample=sample_list)),
        S3.remote(expand(join('data-bucket/gatk_calling', outdir, '02_variants_reference_exome/unfiltered/{sample}_unfiltered.vcf'), sample=sample_list)),
        S3.remote(expand(join('data-bucket/gatk_calling', outdir, '02_variants_reference_exome/filtered/{sample}_filtered.vcf'), sample=sample_list)),
        # S3.remote(expand(join('data-bucket/gatk_calling', outdir, '03_variants_TvN/unfiltered/{patient_tp}_unfiltered.vcf'), patient_tp=tn_pt_tps)),
        # S3.remote(expand(join('data-bucket/gatk_calling', outdir, '03_variants_TvN/filtered/{patient_tp}_filtered.vcf'), patient_tp=tn_pt_tps)),
        S3.remote(expand(join('data-bucket/gatk_calling', outdir, '03_variants_TvN_exome/unfiltered/{patient_tp}_unfiltered.vcf'), patient_tp=tn_pt_tps)),
        S3.remote(expand(join('data-bucket/gatk_calling', outdir, '03_variants_TvN_exome/filtered/{patient_tp}_filtered.vcf'), patient_tp=tn_pt_tps)),


# call mutations with mutect2
# this is just the tumor vs reference analysis
# limit to the segment of regions we have seq data for
rule mutect2_reference:
    input: 
        bam = S3.remote(join('data-bucket/gatk_calling', outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam')),
        bai = S3.remote(join('data-bucket/gatk_calling', outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam.bai')), 
        ref = REF_FILE,
        targets = idt_target_w100
    output: 
        vcf = join(outdir, '02_variants_reference/unfiltered/{sample}_unfiltered.vcf'),
        stats = join(outdir, '02_variants_reference/unfiltered/{sample}_unfiltered.vcf.stats'),
        f1r2 = join(outdir, '02_variants_reference/unfiltered/{sample}_f1r2.tar.gz')
    threads: 8
    shell: """
        gatk Mutect2 \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            -L {input.targets} \
            --native-pair-hmm-threads {threads} \
            --f1r2-tar-gz {output.f1r2}
    """

# filter mutect2 calls
rule mutect2_reference_filter:
    input: 
        vcf = rules.mutect2_reference.output.vcf,
        stats = rules.mutect2_reference.output.stats,
        ref = REF_FILE,
    output: 
        vcf = join(outdir, '02_variants_reference/filtered/{sample}_filtered.vcf'),
    shell: """
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            -stats {input.stats} \
            -O {output.vcf}
    """

# copy target vcfs to remote
rule copy_vcf_remote:
    input: 
        vcf_uf = rules.mutect2_reference.output.vcf,
        stats_uf = rules.mutect2_reference.output.stats,
        f1r2_uf = rules.mutect2_reference.output.f1r2,
        vcf_f = rules.mutect2_reference_filter.output.vcf,
    output:
        vcf_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '02_variants_reference/unfiltered/{sample}_unfiltered.vcf')),
        stats_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '02_variants_reference/unfiltered/{sample}_unfiltered.vcf.stats')),
        f1r2_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '02_variants_reference/unfiltered/{sample}_f1r2.tar.gz')),
        vcf_f = S3.remote(join('data-bucket/gatk_calling', outdir, '02_variants_reference/filtered/{sample}_filtered.vcf')),
    shell: """
        cp {input.vcf_uf} {output.vcf_uf}
        cp {input.stats_uf} {output.stats_uf}
        cp {input.f1r2_uf} {output.f1r2_uf}
        cp {input.vcf_f} {output.vcf_f}
    """

# call mutations with mutect2
# this is just the tumor vs reference analysis
# limit to the segment of regions we have seq data for
rule mutect2_reference_exome:
    input: 
        bam = S3.remote(join('data-bucket/gatk_calling', outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam')),
        bai = S3.remote(join('data-bucket/gatk_calling', outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam.bai')),
        ref = REF_FILE,
    output: 
        vcf = join(outdir, '02_variants_reference_exome/unfiltered/{sample}_unfiltered.vcf'),
        stats = join(outdir, '02_variants_reference_exome/unfiltered/{sample}_unfiltered.vcf.stats'),
        f1r2 = join(outdir, '02_variants_reference_exome/unfiltered/{sample}_f1r2.tar.gz')
    threads: 8
    shell: """
        gatk Mutect2 \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            --native-pair-hmm-threads {threads} \
            --f1r2-tar-gz {output.f1r2}
    """

# filter mutect2 calls
rule mutect2_reference_exome_filter:
    input: 
        vcf = rules.mutect2_reference_exome.output.vcf,
        stats = rules.mutect2_reference_exome.output.stats,
        ref = REF_FILE,
    output: 
        vcf = join(outdir, '02_variants_reference_exome/filtered/{sample}_filtered.vcf'),
    shell: """
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            -stats {input.stats} \
            -O {output.vcf}
    """

# copy target vcfs to remote
rule copy_exome_vcf_remote:
    input: 
        vcf_uf = rules.mutect2_reference_exome.output.vcf,
        stats_uf = rules.mutect2_reference_exome.output.stats,
        f1r2_uf = rules.mutect2_reference_exome.output.f1r2,
        vcf_f = rules.mutect2_reference_exome_filter.output.vcf,
    output:
        vcf_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '02_variants_reference_exome/unfiltered/{sample}_unfiltered.vcf')),
        stats_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '02_variants_reference_exome/unfiltered/{sample}_unfiltered.vcf.stats')),
        f1r2_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '02_variants_reference_exome/unfiltered/{sample}_f1r2.tar.gz')),
        vcf_f = S3.remote(join('data-bucket/gatk_calling', outdir, '02_variants_reference_exome/filtered/{sample}_filtered.vcf')),
    shell: """
        cp {input.vcf_uf} {output.vcf_uf}
        cp {input.stats_uf} {output.stats_uf}
        cp {input.f1r2_uf} {output.f1r2_uf}
        cp {input.vcf_f} {output.vcf_f}
    """

#####################################################################################################
### TUMOR VS NORMAL CALLING #########################################################################
#####################################################################################################
# run this when we have at least one tumor and normal dataset for each patient
# then do the same annotations 
rule mutect2_TvN:
    input: 
        bams = lambda wildcards: S3.remote(expand(join('data-bucket/gatk_calling', outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'), sample=pt_tp_to_samples[wildcards.patient_tp])),
        bais = lambda wildcards: S3.remote(expand(join('data-bucket/gatk_calling', outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam.bai'), sample=pt_tp_to_samples[wildcards.patient_tp])),
        ref = REF_FILE,
        targets = idt_target_w100,
    output:
        vcf = join(outdir, '03_variants_TvN/unfiltered/{patient_tp}_unfiltered.vcf'),
        stats = join(outdir, '03_variants_TvN/unfiltered/{patient_tp}_unfiltered.vcf.stats'),
        f1r2 = join(outdir, '03_variants_TvN/unfiltered/{patient_tp}_f1r2.tar.gz')
    params: 
        sample_bam_string = lambda wildcards: ' -I '.join(pt_tp_to_final_bams[wildcards.patient_tp]),
        normal_bam_string = lambda wildcards: ' -normal '.join(n_samples_pt_tp[wildcards.patient_tp]),
    threads: 8
    shell: """
        gatk Mutect2 \
            -R {input.ref} \
            -I {params.sample_bam_string} \
            -normal {params.normal_bam_string} \
            -O {output.vcf} \
            -L {input.targets} \
            --native-pair-hmm-threads {threads} \
            --f1r2-tar-gz {output.f1r2}
    """


# Read orientation bias detection 
rule mutect2_TvN_readOrientation:
    input: 
        f1r2 = rules.mutect2_TvN.output.f1r2,
    output: 
        model = join(outdir, '03_variants_TvN/unfiltered/{patient_tp}_readOrientation_model.tar.gz'),
    shell: """
        gatk LearnReadOrientationModel -I {input} -O {output}
    """

# filter mutect2 calls
rule mutect2_TvN_filter:
    input: 
        vcf = rules.mutect2_TvN.output.vcf,
        stats = rules.mutect2_TvN.output.stats,
        ref = REF_FILE,
        readOrientation = rules.mutect2_TvN_readOrientation.output,
        # CalculateContamination = lambda wildcards: expand(join(outdir, '01_prepare_bam/CalculateContamination/{t_sample}_contam.table'), t_sample=t_samples_pt_tp[pt_tp])
    output: 
        vcf = join(outdir, '03_variants_TvN/filtered/{patient_tp}_filtered.vcf'),
    params:
        # contamination_string = lambda wildcards: ' --contamination-table '.join(expand(join(outdir, '01_prepare_bam/CalculateContamination/{t_sample}_contam.table'), t_sample=t_samples_pt_tp[pt_tp]))
    shell: """
        # new version with the contaminant info
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            -stats {input.stats} \
            --ob-priors {input.readOrientation} \
            -O {output.vcf}
    """
    # not using contamination calcs
    # --contamination-table {params.contamination_string} \

# copy target vcfs to remote
rule copy_TvN_vcf_remote:
    input: 
        vcf_uf = rules.mutect2_TvN.output.vcf,
        stats_uf = rules.mutect2_TvN.output.stats,
        f1r2_uf = rules.mutect2_TvN.output.f1r2,
        vcf_f = rules.mutect2_TvN_filter.output.vcf,
    output:
        vcf_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '03_variants_TvN/unfiltered/{patient_tp}_unfiltered.vcf')),
        stats_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '03_variants_TvN/unfiltered/{patient_tp}_unfiltered.vcf.stats')),
        f1r2_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '03_variants_TvN/unfiltered/{patient_tp}_f1r2.tar.gz')),
        vcf_f = S3.remote(join('data-bucket/gatk_calling', outdir, '03_variants_TvN/filtered/{patient_tp}_filtered.vcf')),
        model_f = S3.remote(join('data-bucket/gatk_calling', outdir, '03_variants_TvN/unfiltered/{patient_tp}_readOrientation_model.tar.gz')),
    shell: """
        cp {input.vcf_uf} {output.vcf_uf}
        cp {input.stats_uf} {output.stats_uf}
        cp {input.f1r2_uf} {output.f1r2_uf}
        cp {input.vcf_f} {output.vcf_f}
    """

rule mutect2_TvN_exome:
    input: 
        bams = lambda wildcards: S3.remote(expand(join('data-bucket/gatk_calling', outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'), sample=pt_tp_to_samples[wildcards.patient_tp])),
        bais = lambda wildcards: S3.remote(expand(join('data-bucket/gatk_calling', outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam.bai'), sample=pt_tp_to_samples[wildcards.patient_tp])),
        ref = REF_FILE,
    output:
        vcf = join(outdir, '03_variants_TvN_exome/unfiltered/{patient_tp}_unfiltered.vcf'),
        stats = join(outdir, '03_variants_TvN_exome/unfiltered/{patient_tp}_unfiltered.vcf.stats'),
        f1r2 = join(outdir, '03_variants_TvN_exome/unfiltered/{patient_tp}_f1r2.tar.gz')
    params: 
        sample_bam_string = lambda wildcards: ' -I '.join(pt_tp_to_final_bams[wildcards.patient_tp]),
        normal_bam_string = lambda wildcards: ' -normal '.join(n_samples_pt_tp[wildcards.patient_tp]),
    threads: 8
    shell: """
        gatk Mutect2 \
            -R {input.ref} \
            -I {params.sample_bam_string} \
            -normal {params.normal_bam_string} \
            -O {output.vcf} \
            --native-pair-hmm-threads {threads} \
            --f1r2-tar-gz {output.f1r2}
    """


# Read orientation bias detection 
rule mutect2_TvN_exome_readOrientation:
    input: 
        f1r2 = rules.mutect2_TvN_exome.output.f1r2,
    output: 
        model = join(outdir, '03_variants_TvN_exome/unfiltered/{patient_tp}_readOrientation_model.tar.gz'),
    shell: """
        gatk LearnReadOrientationModel -I {input} -O {output}
    """

# filter mutect2 calls
rule mutect2_TvN_exome_filter:
    input: 
        vcf = rules.mutect2_TvN_exome.output.vcf,
        stats = rules.mutect2_TvN_exome.output.stats,
        ref = REF_FILE,
        readOrientation = rules.mutect2_TvN_exome_readOrientation.output,
        # CalculateContamination = lambda wildcards: expand(join(outdir, '01_prepare_bam/CalculateContamination/{t_sample}_contam.table'), t_sample=t_samples_pt_tp[pt_tp])
    output: 
        vcf = join(outdir, '03_variants_TvN_exome/filtered/{patient_tp}_filtered.vcf'),
    params:
        # contamination_string = lambda wildcards: ' --contamination-table '.join(expand(join(outdir, '01_prepare_bam/CalculateContamination/{t_sample}_contam.table'), t_sample=t_samples_pt_tp[pt_tp]))
    shell: """
        # new version with the contaminant info
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            -stats {input.stats} \
            --ob-priors {input.readOrientation} \
            -O {output.vcf}
    """
    # not using contamination calcs
    # --contamination-table {params.contamination_string} \

# copy exome TvN vcfs to remote
rule copy_TvN_exome_vcf_remote:
    input: 
        vcf_uf = rules.mutect2_TvN_exome.output.vcf,
        stats_uf = rules.mutect2_TvN_exome.output.stats,
        f1r2_uf = rules.mutect2_TvN_exome.output.f1r2,
        vcf_f = rules.mutect2_TvN_exome_filter.output.vcf,
    output:
        vcf_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '03_variants_TvN_exome/unfiltered/{patient_tp}_unfiltered.vcf')),
        stats_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '03_variants_TvN_exome/unfiltered/{patient_tp}_unfiltered.vcf.stats')),
        f1r2_uf = S3.remote(join('data-bucket/gatk_calling', outdir, '03_variants_TvN_exome/unfiltered/{patient_tp}_f1r2.tar.gz')),
        vcf_f = S3.remote(join('data-bucket/gatk_calling', outdir, '03_variants_TvN_exome/filtered/{patient_tp}_filtered.vcf')),
        model_f = S3.remote(join('data-bucket/gatk_calling', outdir, '03_variants_TvN_exome/unfiltered/{patient_tp}_readOrientation_model.tar.gz')),
    shell: """
        cp {input.vcf_uf} {output.vcf_uf}
        cp {input.stats_uf} {output.stats_uf}
        cp {input.f1r2_uf} {output.f1r2_uf}
        cp {input.vcf_f} {output.vcf_f}
    """

