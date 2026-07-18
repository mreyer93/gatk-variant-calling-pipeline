#####################################################################################################
### SAMPLE VS REFERENCE CALLING #####################################################################
#####################################################################################################
# use empty string for targets if targets file is not specified
if target_bed_w100 != '':
    target_input_list = [target_bed_w100]
    targets_string = '-L ' + target_bed_w100
else:
    target_input_list = []
    targets_string = ''


# call mutations with mutect2
# this is just the tumor vs reference analysis
# limit to the segment of regions we have seq data for
rule mutect2_reference:
    input:
        target_input_list,
        bam = lambda wildcards: get_final_bam(wildcards.sample),
        bai = lambda wildcards: get_final_bam(wildcards.sample) + '.bai',
        ref = REF_FILE,
        gnomad = config['gnomad_file'],
    output:
        vcf = join(outdir, '02_variants_reference/01_gatk_variant_calling/unfiltered/{sample}_unfiltered.vcf'),
        stats = join(outdir, '02_variants_reference/01_gatk_variant_calling/unfiltered/{sample}_unfiltered.vcf.stats'),
        f1r2 = join(outdir, '02_variants_reference/01_gatk_variant_calling/unfiltered/{sample}_f1r2.tar.gz')
    threads: 4
    shell: """
        gatk Mutect2 \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            {targets_string} \
            --germline-resource {input.gnomad} \
            --native-pair-hmm-threads {threads} \
            --f1r2-tar-gz {output.f1r2}
        # fix error with M chromosome showing up, when it should be chrM
        sed -i "s/contig=<ID=M/contig=<ID=chrM/g" {output.vcf}
    """

# Read orientation bias detection (this sample has no matched normal, so there's no
# contamination-table step here - see variant_calling_TvN.smk for that path)
rule mutect2_reference_readOrientation:
    input:
        f1r2 = rules.mutect2_reference.output.f1r2,
    output:
        model = join(outdir, '02_variants_reference/01_gatk_variant_calling/unfiltered/{sample}_readOrientation_model.tar.gz'),
    shell: """
        gatk LearnReadOrientationModel -I {input} -O {output}
    """

# filter mutect2 calls
rule mutect2_reference_filter:
    input:
        vcf = rules.mutect2_reference.output.vcf,
        stats = rules.mutect2_reference.output.stats,
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
