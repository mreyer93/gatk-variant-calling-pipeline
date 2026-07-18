#####################################################################################################
### TUMOR VS NORMAL CALLING #########################################################################
#####################################################################################################
# use empty string for targets if targets file is not specified
if target_bed_w100 != '':
    target_input_list = [target_bed_w100]
    targets_string = '-L ' + target_bed_w100
else:
    target_input_list = []
    targets_string = ''

# run this when we have at least one tumor and normal dataset for each patient
rule mutect2_TvN:
    input:
        target_input_list,
        bams = lambda wildcards: [get_final_bam(sample) for sample in pt_tp_to_samples[wildcards.sample]],
        ref = REF_FILE,
        gnomad = config['gnomad_file'],
    output:
        vcf = join(outdir, '03_variants_TvN/01_gatk_variant_calling/unfiltered/{sample}_unfiltered.vcf'),
        stats = join(outdir, '03_variants_TvN/01_gatk_variant_calling/unfiltered/{sample}_unfiltered.vcf.stats'),
        f1r2 = join(outdir, '03_variants_TvN/01_gatk_variant_calling/unfiltered/{sample}_f1r2.tar.gz')
    params:
        sample_bam_string = lambda wildcards: ' -I '.join(pt_tp_to_final_bams[wildcards.sample]),
        normal_bam_string = lambda wildcards: ' -normal '.join(n_samples_pt_tp[wildcards.sample]),
    threads: 4
    shell: """
        gatk Mutect2 \
            -R {input.ref} \
            -I {params.sample_bam_string} \
            -normal {params.normal_bam_string} \
            -O {output.vcf} \
            {targets_string} \
            --germline-resource {input.gnomad} \
            --native-pair-hmm-threads {threads} \
            --f1r2-tar-gz {output.f1r2}
        # fix error with M chromosome showing up, when it should be chrM
        sed -i "s/contig=<ID=M/contig=<ID=chrM/g" {output.vcf}

    """

# Read orientation bias detection 
rule mutect2_TvN_readOrientation:
    input: 
        f1r2 = rules.mutect2_TvN.output.f1r2,
    output: 
        model = join(outdir, '03_variants_TvN/01_gatk_variant_calling/unfiltered/{sample}_readOrientation_model.tar.gz'),
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
