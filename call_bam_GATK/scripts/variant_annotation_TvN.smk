#####################################################################################################
### TUMOR VS NORMAL CALLING #########################################################################
#####################################################################################################
# CADD score annotation of the variants against the reference
rule CADD_TvN:
    input:
        join(outdir, '03_variants_TvN/01_gatk_variant_calling/filtered/{patient_tp}_filtered.vcf'),
    output:
        join(outdir, '03_variants_TvN/02_variant_annotations/annotate_CADD/full_table/{patient_tp}_filtered_CADD.vcf.gz')
    params:
        tmpfile = join(outdir, '{patient_tp}_tmp.vcf'),
        genome_version = genome_version
    threads: 1
    shell: """
        sed "s/^chr//g" {input} > {params.tmpfile}
        CADD.sh -a -c {threads} -g {params.genome_version} -o {output} {params.tmpfile}
        rm {params.tmpfile}
    """

rule subset_columns_TvN:
    input:
        vcf = join(outdir, '03_variants_TvN/02_variant_annotations/annotate_CADD/full_table/{patient_tp}_filtered_CADD.vcf.gz')
    output:
        vcf = join(outdir, '03_variants_TvN/02_variant_annotations/annotate_CADD/simple_table/{patient_tp}_filtered_CADD_simple.vcf')
    shell: """
        zcat {input} | cut -f 1,2,3,4,5,6,10,21,116 > {output}
        """

# annotation of variants vs reference with the variantAnnotator
# NO LONGER USED, SUPERCEDED BY FUNCOTATOR
'''
rule variantAnnotator_TvN:
    input: 
        vcf = join(outdir, '03_variants_TvN/01_gatk_variant_calling/filtered/{patient_tp}_filtered.vcf'),
        bam = join(outdir, '01_prepare_bam/recalibrate/{sample}.fixchr.bam'),
        ref = REF_FILE,
        dbsnp = dbsnp_file,
        targets = target_bed_w100,
    output: 
        vcf = join(outdir, '03_variants_TvN/annotate_VariantAnnotator/{patient_tp}_filtered_annot.vcf'),
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
rule funcotator_TvN:
    input: 
        vcf = join(outdir, '03_variants_TvN/01_gatk_variant_calling/filtered/{patient_tp}_filtered.vcf'),
        ref = REF_FILE,
        dbsnp = dbsnp_file,
    output: 
        vcf = join(outdir, '03_variants_TvN/02_variant_annotations/annotate_Funcotator/{patient_tp}_filtered_fucnc.vcf'),
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
rule combine_annoations_TvN:
    input:
        vcf_funcotator = rules.funcotator_TvN.output.vcf, 
        vcf_cadd = rules.subset_columns_TvN.output.vcf,
    output: 
        vcf_unfilt = join(outdir, '03_variants_TvN/02_variant_annotations/annotations_combined/unfiltered/{patient_tp}_unfiltered_annotated.vcf'),
        vcf_filt = join(outdir, '03_variants_TvN/02_variant_annotations/annotations_combined/filtered/{patient_tp}_filtered_annotated.vcf')
    params: 
        filter_min_af = config['filter_min_af'],
        filter_min_dp = config['filter_min_dp'],
        filter_GATK = config['filter_GATK'],
        filter_min_cadd = config['filter_min_cadd']
    script: "combine_annotation_TvN.R"
