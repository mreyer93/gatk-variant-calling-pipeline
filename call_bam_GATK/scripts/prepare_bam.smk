# steps to prepare the bam files for variant calling

# first, add read groups to the samples that need it
rule add_read_groups:
    input: lambda wildcards: bam_map[wildcards.sample]
    output: temporary(join(outdir, '01_prepare_bam/add_readgroups/{sample}.bam'))
    params:
        RGSM = lambda wildcards: wildcards.sample,
        RGLB = lambda wildcards: metadata.loc[wildcards.sample, 'library'],
        RGPL = lambda wildcards: metadata.loc[wildcards.sample, 'platform'],
        RGPU = lambda wildcards: metadata.loc[wildcards.sample, 'platform_unit'],
    shell: """
        picard AddOrReplaceReadGroups \
            -I {input} \
            -O {output} \
            --RGSM {params.RGSM} \
            --RGLB {params.RGLB} \
            --RGPL {params.RGPL} \
            --RGPU {params.RGPU}
    """

# fix chromosome annotations: adds chr to all chromosomes in the bam files so 
# the annotation pipelines work correctly
rule fix_chromosomes:
    input: rules.add_read_groups.output
    output: 
        bam = temporary(join(outdir, '01_prepare_bam/add_readgroups/{sample}.fixchr.bam')),
        bai = temporary(join(outdir, '01_prepare_bam/add_readgroups/{sample}.fixchr.bam.bai'))
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
            sed -e 's/SN:MT/SN:chrM/' | samtools reheader - {input} > {output.bam}
        samtools index {output.bam}
    """

# base quality score recalibration
rule bqsr:
    input: 
        bam = rules.fix_chromosomes.output.bam,
        bai = rules.fix_chromosomes.output.bai,
        ref = REF_FILE,
        dbsnp = dbsnp_file,
        targets = target_bed_w100
    output: 
        table = join(outdir, '01_prepare_bam/recalibrate/{sample}.table'),
        bam = join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'),
        bai = join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam.bai')
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
        samtools index {output.bam}
    """

# Calculate depth based on targets 
rule depth:
    input: 
        bam = lambda wildcards: get_final_bam(wildcards.sample),
        targets = targets_depth,
    output: 
        join(outdir, '01_prepare_bam/depth/{sample}.tsv.gz'),
    shell: """
        samtools depth \
        -d 0 \
        -a \
        -b {input.targets} \
        {input.bam} | gzip > {output}
        """

# smooth depth result
rule smooth_depth:
    input: 
        d = rules.depth.output,
        targets = targets_depth,
    output: 
        join(outdir, '01_prepare_bam/depth_smoothed/{sample}.tsv'),
    script: "smooth_depth.R"

# join smooth depth result
rule join_smooth_depth:
    input: 
        expand(join(outdir, '01_prepare_bam/depth_smoothed/{sample}.tsv'), sample=sample_list)
    output: 
        join(outdir, '01_prepare_bam/depth_smoothed_join.tsv')
    params:
        sample_list = sample_list
    script: "join_smoothed_depth.R"

# Run GetPileupSummaries to summarize read support for a set number of known variant sites.
# gets run on EVERY Bam
rule GetPileupSummaries:
    input: 
        bam = lambda wildcards: get_final_bam(wildcards.sample),
        targets = target_bed_w100,
        gnomad = config['gnomad_file'],
    output: 
        table = join(outdir, '01_prepare_bam/GetPileupSummaries/{sample}_GetPileupSummaries.table'),
    shell: """
        gatk GetPileupSummaries \
            -I {input.bam} \
            -V {input.gnomad} \
            -L {input.targets} \
            -O {output}
        """

# Run CalculateContamination to understand tumor normal contamination
# run this for each tumor, normal pair
rule CalculateContamination:
    input: 
        tumor_table = lambda wildcards: join(outdir, '01_prepare_bam/GetPileupSummaries/' + pt_contamination_tn_sets[wildcards.t_sample][0] + '_GetPileupSummaries.table'),
        normal_table = lambda wildcards: join(outdir, '01_prepare_bam/GetPileupSummaries/' + pt_contamination_tn_sets[wildcards.t_sample][1] + '_GetPileupSummaries.table'),
    output: 
        table = join(outdir, '01_prepare_bam/CalculateContamination/{t_sample}_contam.table'),
    shell: """
        gatk CalculateContamination \
           -I {input.tumor_table} \
           -matched {input.normal_table} \
           -O {output}
        """
