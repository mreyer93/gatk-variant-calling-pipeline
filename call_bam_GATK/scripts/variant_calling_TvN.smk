#####################################################################################################
### TUMOR VS NORMAL CALLING #########################################################################
#####################################################################################################
# Mutect2 runs once per interval shard (scripts/scatter.smk), limited to the targets the
# shard covers; the shard VCFs, stats and F1R2 counts are gathered below.
m2_tvn_dir = join(outdir, '03_variants_TvN/01_gatk_variant_calling/unfiltered')

# run this when we have at least one tumor and normal dataset for each patient
rule mutect2_TvN:
    input:
        bams = lambda wildcards: [get_final_bam(sample) for sample in pt_tp_to_samples[wildcards.sample]],
        ref = REF_FILE,
        gnomad = config['gnomad_file'],
        intervals = join(outdir, '00_scatter/{shard}.bed'),
    output:
        vcf = temp(join(m2_tvn_dir, 'scatter/{shard}/{sample}.vcf')),
        stats = temp(join(m2_tvn_dir, 'scatter/{shard}/{sample}.vcf.stats')),
        f1r2 = temp(join(m2_tvn_dir, 'scatter/{shard}/{sample}_f1r2.tar.gz')),
    params:
        sample_bam_string = lambda wildcards: ' -I '.join(pt_tp_to_final_bams[wildcards.sample]),
        normal_bam_string = lambda wildcards: ' -normal '.join(n_samples_pt_tp[wildcards.sample]),
    threads: gatk_threads
    resources: mem_mb = gatk_java_mem_gb * 1024 + 1024
    shell: """
        gatk --java-options "-Xmx{gatk_java_mem_gb}g" Mutect2 \
            -R {input.ref} \
            -I {params.sample_bam_string} \
            -normal {params.normal_bam_string} \
            -O {output.vcf} \
            -L {input.intervals} \
            --germline-resource {input.gnomad} \
            --native-pair-hmm-threads {threads} \
            --f1r2-tar-gz {output.f1r2}
    """

# gather the shards: calls, and the callable-site stats FilterMutectCalls needs
rule mutect2_TvN_gather:
    input:
        vcfs = expand(join(m2_tvn_dir, 'scatter/{shard}/{{sample}}.vcf'), shard=SHARDS),
        stats = expand(join(m2_tvn_dir, 'scatter/{shard}/{{sample}}.vcf.stats'), shard=SHARDS),
    output:
        vcf = join(m2_tvn_dir, '{sample}_unfiltered.vcf'),
        stats = join(m2_tvn_dir, '{sample}_unfiltered.vcf.stats'),
    params:
        vcfs = lambda wildcards, input: ' '.join('-I ' + f for f in input.vcfs),
        stats = lambda wildcards, input: ' '.join('-stats ' + f for f in input.stats),
    shell: """
        gatk MergeVcfs {params.vcfs} -O {output.vcf}
        gatk MergeMutectStats {params.stats} -O {output.stats}
        # fix error with M chromosome showing up, when it should be chrM.
        # Written via a temp file rather than `sed -i`: BSD/macOS sed requires an
        # explicit backup suffix after -i, so `sed -i "expr" file` fails there.
        sed "s/contig=<ID=M/contig=<ID=chrM/g" {output.vcf} > {output.vcf}.tmp \
            && mv {output.vcf}.tmp {output.vcf}
    """

# Read orientation bias detection, from every shard's F1R2 counts together
rule mutect2_TvN_readOrientation:
    input:
        f1r2 = expand(join(m2_tvn_dir, 'scatter/{shard}/{{sample}}_f1r2.tar.gz'), shard=SHARDS),
    output:
        model = join(m2_tvn_dir, '{sample}_readOrientation_model.tar.gz'),
    params:
        f1r2 = lambda wildcards, input: ' '.join('-I ' + f for f in input.f1r2),
    shell: """
        gatk LearnReadOrientationModel {params.f1r2} -O {output}
    """

# filter mutect2 calls
rule mutect2_TvN_filter:
    input:
        vcf = rules.mutect2_TvN_gather.output.vcf,
        stats = rules.mutect2_TvN_gather.output.stats,
        ref = REF_FILE,
        readOrientation = rules.mutect2_TvN_readOrientation.output,
        CalculateContamination = lambda wildcards: expand(join(outdir, '01_prepare_bam/CalculateContamination/{t_sample}_contam.table'), t_sample=t_samples_pt_tp[wildcards.sample])
    output:
        vcf = join(outdir, '03_variants_TvN/01_gatk_variant_calling/filtered/{sample}_filtered.vcf'),
    params:
        contamination_string = lambda wildcards: ' --contamination-table '.join(expand(join(outdir, '01_prepare_bam/CalculateContamination/{t_sample}_contam.table'), t_sample=t_samples_pt_tp[wildcards.sample]))
    shell: """
        # new version with the contaminant info
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            -stats {input.stats} \
            --contamination-table {params.contamination_string} \
            --ob-priors {input.readOrientation} \
            -O {output.vcf}
    """
