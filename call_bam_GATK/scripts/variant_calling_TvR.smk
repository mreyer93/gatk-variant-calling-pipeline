#####################################################################################################
### SAMPLE VS REFERENCE CALLING #####################################################################
#####################################################################################################
# Mutect2 runs once per interval shard (scripts/scatter.smk), limited to the targets the
# shard covers; the shard VCFs, stats and F1R2 counts are gathered below.
m2_ref_dir = join(outdir, '02_variants_reference/01_gatk_variant_calling/unfiltered')

# call mutations with mutect2
# this is just the tumor vs reference analysis
# limit to the segment of regions we have seq data for
rule mutect2_reference:
    input:
        bam = lambda wildcards: get_final_bam(wildcards.sample),
        bai = lambda wildcards: get_final_bam(wildcards.sample) + '.bai',
        ref = REF_FILE,
        gnomad = config['gnomad_file'],
        intervals = join(outdir, '00_scatter/{shard}.bed'),
    output:
        vcf = temp(join(m2_ref_dir, 'scatter/{shard}/{sample}.vcf')),
        stats = temp(join(m2_ref_dir, 'scatter/{shard}/{sample}.vcf.stats')),
        f1r2 = temp(join(m2_ref_dir, 'scatter/{shard}/{sample}_f1r2.tar.gz')),
    threads: gatk_threads
    resources: mem_mb = gatk_java_mem_gb * 1024 + 1024
    shell: """
        gatk --java-options "-Xmx{gatk_java_mem_gb}g" Mutect2 \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            -L {input.intervals} \
            --germline-resource {input.gnomad} \
            --native-pair-hmm-threads {threads} \
            --f1r2-tar-gz {output.f1r2}
    """

# gather the shards: calls, and the callable-site stats FilterMutectCalls needs
rule mutect2_reference_gather:
    input:
        vcfs = expand(join(m2_ref_dir, 'scatter/{shard}/{{sample}}.vcf'), shard=SHARDS),
        stats = expand(join(m2_ref_dir, 'scatter/{shard}/{{sample}}.vcf.stats'), shard=SHARDS),
    output:
        vcf = join(m2_ref_dir, '{sample}_unfiltered.vcf'),
        stats = join(m2_ref_dir, '{sample}_unfiltered.vcf.stats'),
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

# Read orientation bias detection, from every shard's F1R2 counts together (this sample
# has no matched normal, so there's no contamination-table step here - see
# variant_calling_TvN.smk for that path)
rule mutect2_reference_readOrientation:
    input:
        f1r2 = expand(join(m2_ref_dir, 'scatter/{shard}/{{sample}}_f1r2.tar.gz'), shard=SHARDS),
    output:
        model = join(m2_ref_dir, '{sample}_readOrientation_model.tar.gz'),
    params:
        f1r2 = lambda wildcards, input: ' '.join('-I ' + f for f in input.f1r2),
    shell: """
        gatk LearnReadOrientationModel {params.f1r2} -O {output}
    """

# filter mutect2 calls
rule mutect2_reference_filter:
    input:
        vcf = rules.mutect2_reference_gather.output.vcf,
        stats = rules.mutect2_reference_gather.output.stats,
        ref = REF_FILE,
        readOrientation = rules.mutect2_reference_readOrientation.output,
    output:
        vcf = join(outdir, '02_variants_reference/01_gatk_variant_calling/filtered/{sample}_filtered.vcf'),
    shell: """
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            -stats {input.stats} \
            --ob-priors {input.readOrientation} \
            -O {output.vcf}
    """
