import pandas as pd
from os.path import join, splitext
from collections import defaultdict

# Setup for the GATK variant calling pipeline

# define parameters from the configfile 
outdir = config['output_directory']
# for now, used the ref file with fixed chromosome names
REF_FILE = config['REF_FILE']
dbsnp_file =config['dbsnp_file']
dbsnp_common_file =config['dbsnp_common_file']
# targets with 100bp windows for mutation calling
if 'target_bed_w100' in config:
    target_bed_w100 = config['target_bed_w100']
else:
    target_bed_w100 = config['targets']
# targets with 100bp windows for depth calling
if 'targets_depth' in config:
    targets_depth = config['targets_depth']
else:
    targets_depth = target_bed_w100
# genome version for CADD scoring
genome_version = config['genome_version'] 

metadata = pd.read_csv(config["bam_metadata"], delimiter = '\t',index_col=False)
req_columns = ['sample', 'bamfile', 'patient', 'timepoint', 'tumor_normal', 'library', 'platform', 'platform_unit']
req_columns_small = ['sample', 'bamfile']
if not(all([r in list(metadata.columns) for r in req_columns])):
    if all([r in list(metadata.columns) for r in req_columns_small]):
        print("###################################################################################")
        print("WARNING: not all expected columns are present, assuming all samples are independent")
        print("###################################################################################")
        metadata['patient'] = list(range(metadata.shape[0]))
        metadata['timepoint'] = 1
        metadata['tumor_normal'] = "T"
        metadata['library'] = "lib1"
        metadata['platform'] = "ILLUMINA"
        metadata['platform_unit'] = "unit1"
    else:
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

# add annotation for patient-timepoint pair
metadata['pt_tp'] = [str(a)+'_'+str(b) for a, b in zip(metadata['patient'], metadata['timepoint'])]       
# sample name to final bamfile
sample_to_final_bam = {a:b for a, b in zip(list(metadata['sample']), expand(join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'), sample= list(metadata['sample'])))}
# map out what we need for calling tumor vs normal for each patient-timepoint pair
# this will be tumors from this timepoint, and the DX normals
tn_pt_tps = []
# map from patient-timepoint to samples and final bamfiles
pt_to_samples = {}
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
    pt_to_samples[p] = list(m_pt['sample'])
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
            pt_tp_to_final_bams[pt_tp] = expand(join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'), sample= pt_tp_to_samples[pt_tp])

# for tumor-normal contamination estimation, which pairs do we want to calculate
# it's each tumor vs one normal, ideally the best one
# for now just use the first? 
pt_contamination_normals = {a:b[0] for a,b in n_samples_pt.items()}
# all normal samples with a matched tumor 
all_normal_samples_with_tumor = [item for sublist in list(n_samples_pt.values()) for item in sublist] 

# this gets ruin for every T sample from each patient 
pt_contamination_tn_sets = {}
for pt,samples in t_samples_pt.items():
    for s in samples:
        pt_contamination_tn_sets[s] = [s, pt_contamination_normals[pt]]

all_t_contamination_samples = [item for sublist in list(t_samples_pt.values()) for item in sublist] 

if not('skip_bqsr' in config):
    config['skip_bqsr'] = False
if not('skip_annotation' in config):
    config['skip_annotation'] = False
# determine final bam for variant calling based on config 
def get_final_bam(sample):
    if config['skip_bqsr']:
        final_bam = join(outdir, '01_prepare_bam/add_readgroups/{sample}.fixchr.bam')
    else:
        final_bam = join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam')
    return(final_bam)