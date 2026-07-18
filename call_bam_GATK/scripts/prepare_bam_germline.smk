# steps to prepare the bam files for variant calling

# first, add read groups to the samples that need it
rule add_read_groups:
    input: lambda wildcards: bam_map[wildcards.sample]
    output: temporary(join(outdir, '01_prepare_bam/add_readgroups/{sample}.bam'))
    params:
        RGSM = lambda wildcards: wildcards.sample
    shell: """
        picard AddOrReplaceReadGroups \
            -I {input} \
            -O {output} \
            --RGSM {params.RGSM} \
            --RGLB l1 \
            --RGPL pl1 \
            --RGPU pu1
    """

# fix chromosome annotations: adds chr to all chromosomes in the bam files so 
# the annotation pipelines work correctly
rule fix_chromosomes:
    input: rules.add_read_groups.output
    output: temporary(join(outdir, '01_prepare_bam/add_readgroups/{sample}.fixchr.bam'))
    shell: """
        samtools view -H {input} |\
            sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' | \
            sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' | \
            sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' | \
            sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' | \
            sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | \
            sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' | \
            sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' | \
            sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' | \
            sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' | \
            sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' | \
            sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' | \
            sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' | \
            sed -e 's/SN:MT/SN:chrM/' | samtools reheader - {input} > {output}
    """

rule bam_index:
    input: rules.fix_chromosomes.output
    output: temporary(join(outdir, '01_prepare_bam/add_readgroups/{sample}.fixchr.bam.bai'))
    shell: """
        samtools index {input}
    """

# base quality score recalibration
rule bqsr:
    input: 
        bam = rules.fix_chromosomes.output,
        bai = rules.bam_index.output,
        ref = REF_FILE,
        dbsnp = dbsnp_file,
        targets = targets
    output: 
        table = join(outdir, '01_prepare_bam/recalibrate/{sample}.table'),
        bam = join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'),
        bai = join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bai')
    shell: """
        gatk BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            --known-sites {input.dbsnp} \
            -O {output.table} \
            -L {input.targets} 

        gatk ApplyBQSR \
            -R {input.ref} \
            -I {input.bam} \
            --bqsr-recal-file {output.table} \
            -O {output.bam}
    """

