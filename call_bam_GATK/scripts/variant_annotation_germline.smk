#####################################################################################################
### GERMLINE CALLING ANNOTATION #####################################################################
#####################################################################################################
rule CADD:
    input: 
        rules.hardFilter_select.output,
    output:
        join(outdir, '07_joint_vcf/02_variant_annotations/CADD_full_table.vcf.gz')
    params:
        tmpsplit = join(outdir, 'tmp_split'),
        tmpout = join(outdir, 'tmp_cadd')
    threads: 16
    shell: """
        # split to multiprocess
        # TODO: do this with snakemake rules instead of an xargs so we can resume
        rm -rf {params.tmpsplit}
        rm -rf {params.tmpout}
        mkdir {params.tmpsplit}
        mkdir {params.tmpout}
        
        if ! command -v CADD.sh &> /dev/null; then COMMAND=cadd.sh; else COMMAND=CADD.sh; fi

        echo "SPLITTING!"
        cut -f 1-9 {input} | sed "s/^chr//g" | split --lines 10000 --additional-suffix=.vcf - {params.tmpsplit}/x

        cd {params.tmpsplit}
        ls x*.vcf | xargs -P {threads} -I {{}} sh -c "$COMMAND -a -c 1 -g GRCh37 -o {params.tmpout}/{{}} {{}}.gz"

        zcat *.vcf.gz | pigz -p {threads} > {output}
    """

rule subset_columns:
    input: rules.CADD.output
    output:
        vcf = join(outdir, '07_joint_vcf/02_variant_annotations/CADD_simple_table.vcf')
    shell: """
        zcat {input} | cut -f 1,2,3,4,5,6,10,21,116 > {output}
        """

# annotation of variants vs reference with the funcotator
rule funcotator:
    input: 
        vcf = rules.hardFilter_select.output,
        ref = REF_FILE,
    output: 
        vcf = join(outdir, '07_joint_vcf/02_variant_annotations/funcotator.vcf')
    params:
        funcotator_data_path = config['funcotator_data_path']
    threads: 4
    shell: """
        gatk Funcotator \
            --variant {input.vcf} \
            --reference {input.ref} \
            --ref-version hg19 \
            --data-sources-path {params.funcotator_data_path}  \
            --output {output.vcf} \
            --output-file-format VCF
        """

# combine the annotation from CADD and funcotator
rule combine_annoations:
    input:
        vcf_funcotator = rules.funcotator.output.vcf, 
        vcf_cadd = rules.subset_columns.output.vcf,
    output: 
        vcf = join(outdir, '07_joint_vcf/02_variant_annotations/annotations_combined.vcf'),
    params: 
        filter_min_dp = config['filter_min_dp'],
    script: "combine_annotation_germline.R"
