'''
Sequence data processing pipeline

This pipeline takes human sequencing data and does the following:
 - Quality control and sequence trimming
 - Alignment against the human reference genome
'''

################################################################################
#  import other scripts and functions 
# this first script does configuration and reads the sample file
include: "scripts/wgs_functions.smk"
# creates reference indices if they don't exist
include: "scripts/reference_index.smk"
# fastqc, trimming, multiqc
include: "scripts/qc.smk"
# alignment with bwa and coverage qc
include: "scripts/alignment.smk"

print("     ##################################       ")
print("                   SAMPLES                    ")
print(samples)
print("     ##################################       ")
wildcard_constraints:
    sample='|'.join([x for x in samples]),
################################################################################

rule all:
    input:
        expand(join(outdir,  "00_qc_reports/{step}_multiqc/multiqc_report.html"), step=['pre', 'post']),
        expand(join(outdir, "01_trimmed/{sample}_R1_val_1.fq.gz"), sample=samples),
        expand(join(outdir, "02_align/bam/{sample}.mkdup.bam"), sample=samples),
        join(outdir, '02_align/aligned_counts.txt'),
        expand(join(outdir, "02_align/cov/{sample}.mkdup.bam.w0.cov"), sample=samples),
        expand(join(outdir, "02_align/cov/{sample}.mkdup.bam.w100.cov"), sample=samples),
        expand(join(outdir, "02_align/offtarget/{sample}_off.bam"), sample=samples),
        join(outdir, "02_align/flagstat_offtarget.txt"),
        join(outdir, "primer_check.pdf"),

