#####################################################################################################
### SAMPLE VS REFERENCE CALLING ANNOTATION ##########################################################
#####################################################################################################

# CADD score annotation of the variants against the reference
rule CADD_reference:
    input: join(outdir, '02_variants_reference/01_gatk_variant_calling/filtered/{sample}_filtered.vcf')
    output:
        join(outdir, '02_variants_reference/02_variant_annotations/annotate_CADD/full_table/{sample}_filtered_CADD.vcf.gz')
    params:
        tmpfile = join(outdir, '{sample}_tmp.vcf')
    threads: 1
    shell: """
        sed "s/^chr//g" {input} > {params.tmpfile}
        if ! command -v CADD.sh &> /dev/null
        then
            cadd.sh -a -c {threads} -g {genome_version} -o {output} {params.tmpfile}
        else
            CADD.sh -a -c {threads} -g {genome_version} -o {output} {params.tmpfile}
        fi

        rm {params.tmpfile}
    """

rule subset_columns_reference:
    input: rules.CADD_reference.output
    output:
        vcf = join(outdir, '02_variants_reference/02_variant_annotations/annotate_CADD/simple_table/{sample}_filtered_CADD_simple.vcf')
    shell: """
        zcat {input} | cut -f 1,2,3,4,5,6,10,21,116 > {output}
        """

# annotation of variants vs reference with the variantAnnotator
# NO LONGER USED, SUPERCEDED BY FUNCOTATOR
'''
rule variantAnnotator_reference:
    input: 
        vcf = join(outdir, '02_variants_reference/filtered/{sample}_filtered.vcf'),
        bam = join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'),
        ref = REF_FILE,
        dbsnp = dbsnp_file,
        targets = target_bed_w100,
    output: 
        vcf = join(outdir, '02_variants_reference/annotate_VariantAnnotator/{sample}_filtered_annot.vcf'),
    shell: """
        gatk VariantAnnotator \
            -R {input.ref} \
            -I {input.bam} \
            -V {input.vcf} \
            -O {output.vcf} \
            -A Coverage \
            --dbsnp {input.dbsnp} \
            -L {input.targets} \
        """
'''

# annotation of variants vs reference with the funcotator
rule funcotator_reference:
    input: 
        vcf = join(outdir, '02_variants_reference/01_gatk_variant_calling/filtered/{sample}_filtered.vcf'),
        bam = lambda wildcards: get_final_bam(wildcards.sample),
        ref = REF_FILE,
        dbsnp = dbsnp_file,
    output: 
        vcf = join(outdir, '02_variants_reference/02_variant_annotations/annotate_Funcotator/{sample}_filtered_func.vcf'),
    params:
        funcotator_data_path = config['funcotator_data_path']
    threads: 4
    shell: """
        gatk Funcotator \
            --variant {input.vcf} \
            --reference {input.ref} \
            --ref-version {funcotator_ref_version} \
            --data-sources-path {params.funcotator_data_path}  \
            --output {output.vcf} \
            --output-file-format VCF
        """

# combine the annotation from CADD and funcotator
rule combine_annoations_reference:
    input:
        vcf_funcotator = rules.funcotator_reference.output.vcf, 
        vcf_cadd = rules.subset_columns_reference.output.vcf,
    output: 
        vcf_unfilt = join(outdir, '02_variants_reference/02_variant_annotations/annotations_combined/unfiltered/{sample}_unfiltered_annotated.vcf'),
        vcf_filt = join(outdir, '02_variants_reference/02_variant_annotations/annotations_combined/filtered/{sample}_filtered_annotated.vcf')
    params: 
        filter_min_af = config['filter_min_af'],
        filter_min_dp = config['filter_min_dp'],
        filter_GATK = config['filter_GATK']
    script: "combine_annotation.R"

# combine the annotated VCFs from all samples from each patiet
# in both filtered and unfiltered case
rule combine_sample_vcfs:
    input:
        vcf_unfilt = expand(join(outdir, '02_variants_reference/02_variant_annotations/annotations_combined/unfiltered/{sample}_unfiltered_annotated.vcf'), sample=sample_list),
        vcf_filt = expand(join(outdir, '02_variants_reference/02_variant_annotations/annotations_combined/filtered/{sample}_filtered_annotated.vcf'), sample=sample_list),
        metadata = config["bam_metadata"],
    output: 
        vcf_unfilt = expand(join(outdir, '02_variants_reference/03_variant_new_scores/prescore/unfiltered/{pt}.vcf'), pt=unique_pts),
        vcf_filt = expand(join(outdir, '02_variants_reference/03_variant_new_scores/prescore/filtered/{pt}.vcf'), pt=unique_pts),
        normal_missing_bed = expand(join(outdir, '02_variants_reference/03_variant_new_scores/prescore/to_compute_bedfiles/{sample}.bed'), sample=all_normal_samples_with_tumor)
    params: 
        annot_combined_dir = join(outdir, '02_variants_reference/02_variant_annotations/annotations_combined'),
        outdir_base = join(outdir, '02_variants_reference/03_variant_new_scores/prescore/'),
        filter_min_cadd = config['filter_min_cadd'],
        filter_min_af = config['filter_min_af'],
        filter_min_dp = config['filter_min_dp'],
        filter_GATK = config['filter_GATK']
    script: "combine_samples.R"

# genotype normal samples at the missing sites
rule genotype_missing_normals:
    input:
        normal_missing_bed = join(outdir, '02_variants_reference/03_variant_new_scores/prescore/to_compute_bedfiles/{sample}.bed'),
        bam = lambda wildcards: get_final_bam(wildcards.sample),
        ref = REF_FILE,
    output:
        join(outdir,'02_variants_reference/03_variant_new_scores/prescore/computed_vcfs/{sample}.vcf')
    shell: """
        if [ -s {input.normal_missing_bed} ]; then
            bcftools mpileup -d 99999 -Ov -R {input.normal_missing_bed} -f {input.ref} {input.bam} | bcftools call -c - > {output}
        else
            touch {output}
        fi
    """

rule redo_new_scores:
    input:
        new_scores_filtt = expand(join(outdir, '02_variants_reference/03_variant_new_scores/prescore/filtered/{pt}.vcf'), pt=unique_pts), 
        normal_redone = expand(join(outdir,'02_variants_reference/03_variant_new_scores/prescore/computed_vcfs/{sample}.vcf'), sample=all_normal_samples_with_tumor),
        metadata = config["bam_metadata"],
    output: 
        vcf_redone_filt = expand(join(outdir, '02_variants_reference/03_variant_new_scores/final/filtered/{pt}.vcf'), pt=unique_pts),
        vcf_redone_gatk = expand(join(outdir, '02_variants_reference/03_variant_new_scores/final/filtered_GATK/{pt}.vcf'), pt=unique_pts),
    params: 
        indir_base = join(outdir, '02_variants_reference/03_variant_new_scores/prescore/filtered'),
        indir_vcf = join(outdir, '02_variants_reference/03_variant_new_scores/prescore/computed_vcfs'),
        outdir_base = join(outdir, '02_variants_reference/03_variant_new_scores/final'),
        filter_min_cadd = config['filter_min_cadd'],
        filter_min_af = config['filter_min_af'],
        filter_min_dp = config['filter_min_dp'],
        filter_GATK = config['filter_GATK']
    script: "combine_samples_redo.R"