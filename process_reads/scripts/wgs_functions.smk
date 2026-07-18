###############################################################################
# load python modules
import fnmatch as fn
import pandas as pd
from os.path import join, splitext
from collections import defaultdict
from snakemake.io import Wildcards, expand
import snakemake.rules
import sys

################################################################################
# define parameters from the configfile 
REF_FILE = config['reference_file']
outdir = config['output_directory']
sample_reads_f = config['sample_file']
read_length = config['read_length']
target_bed = config['target_bed']
target_bed_w100 = config['target_bed_w100']
target_bed_w1000 = config['target_bed_w1000']
################################################################################
# read sample read map
sample_reads = pd.read_csv(config["sample_file"], delimiter = '\t', names=['sample', 'read_1', 'read_2'])

samples = list(sample_reads['sample'])
read1_list = list(sample_reads['read_1'])
read2_list = list(sample_reads['read_2'])

# ensure no duplicates in sample, read1, read2
if (any(samples.count(x) > 1  for x in samples)):
    sys.exit("Check sample_reads, duplictes in sample column")
if (any(read1_list.count(x) > 1  for x in read1_list)):
    sys.exit("Check sample_reads, duplictes in read1 column")
if (any(read2_list.count(x) > 1  for x in read2_list)):
    sys.exit("Check sample_reads, duplictes in read2 column")

# ensure all files specified actually exist
check_list = [not os.path.exists(a) for a in read1_list+read2_list]
if(any(check_list)):
    missing_files = [b for a,b in zip(check_list, read1_list + read2_list) if a]
    print("Missing files: ")
    for a in missing_files: print(a) 
    sys.exit("Some specified fastq files do not exist")

# map from sample prefix to reads (list of two)
read_map = {a: [b,c] for a,b,c in zip(samples, read1_list, read2_list)}

# Ensure all files are gzipped
if not all([a.endswith(".gz") for a in read1_list]) and all([a.endswith(".gz") for a in read2_list]):
    sys.exit('All input read files must be gzipped!')
