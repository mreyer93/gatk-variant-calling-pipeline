import pandas as pd
import numpy as np
from os.path import join, splitext
from collections import defaultdict

# Setup for the GATK variant calling pipeline
# list of chromosomes to process 
chromosome_list = [f"chr{i}" for i in list(range(1,23))]

# define parameters from the configfile 
outdir = config['output_directory']
# for now, used the ref file with fixed chromosome names
REF_FILE = config['REF_FILE']
dbsnp_file =config['dbsnp_file']
dbsnp_common_file =config['dbsnp_common_file']
# genome version for CADD scoring (CADD uses GRCh37/GRCh38 naming directly)
genome_version = config['genome_version']
# Funcotator uses UCSC-style hg19/hg38 naming instead of GRCh37/GRCh38
funcotator_ref_version = {'GRCh37': 'hg19', 'GRCh38': 'hg38'}[genome_version]
# targets with 100bp windows for mutation calling
if 'target_bed_w100' in config:
    targets = config['target_bed_w100']
else:
    targets = config['targets']
metadata = pd.read_csv(config["bam_metadata"], delimiter = '\t', index_col=False)
req_columns = ['sample', 'bamfile']
if not(all([r in list(metadata.columns) for r in req_columns])):
    sys.exit("Metadata must contain the following columns " + str(req_columns))
# set rownames to sample
metadata.index = list(metadata['sample'])

# ensure no duplicates in sample, bamfile
sample_list = list(metadata['sample'])
bamfile_list = list(metadata['bamfile'])
if (any(sample_list.count(x) > 1  for x in sample_list)):
    sys.exit("Check metadata, duplictes in sample column")
if (any(bamfile_list.count(x) > 1  for x in bamfile_list)):
    sys.exit("Check metadata, duplictes in bamfile column")
# print(sample_list)
# ensure all files specified actually exist
check_list = [not os.path.exists(a) for a in bamfile_list]
if(any(check_list)) and not config['annotation_only']:
    missing_files = [b for a,b in zip(check_list, bamfile_list) if a]
    print("Missing files: ")
    for a in missing_files: print(a) 
    sys.exit("Some specified bam files do not exist")

# map from sample names to bamfiles
bam_map = {a: b for a,b in zip(sample_list, bamfile_list)}

if not('skip_bqsr' in config):
    config['skip_bqsr'] = False
if not('skip_annotation' in config):
    config['skip_annotation'] = False
if not('final_bam' in config):
    config['final_bam'] = False
# skip_cadd: skip the CADD annotation step (its 400G+ database isn't always available,
# e.g. on smaller local machines). AlphaMissense is a much lighter-weight (~640M) germline
# missense deleteriousness score that can run instead if alphamissense_file is set - see
# manual/germline.md. Unlike CADD it only covers missense SNVs, not indels/nonsense/splice/
# non-coding variants, so those will have no deleteriousness annotation in that mode.
if not('skip_cadd' in config):
    config['skip_cadd'] = False
if not('alphamissense_file' in config):
    config['alphamissense_file'] = ''
# determine final bam for variant calling based on config 
# See the note in setup.smk: the sample name must be substituted here rather than left
# as a literal "{sample}" for Snakemake to re-expand against the requesting rule.
def get_final_bam(sample):
    if config['skip_bqsr']:
        return join(outdir, '01_prepare_bam/add_readgroups/{}.fixchr.bam'.format(sample))
    elif config['final_bam']:
        return bam_map[sample]
    return join(outdir, '01_prepare_bam/recalibrate/{}.fixchr.bam'.format(sample))