#####################################################################################################
### GERMLINE CALLING ANNOTATION #####################################################################
#####################################################################################################
rule alphamissense:
    input:
        vcf = rules.hardFilter_select.output,
        am_file = config['alphamissense_file'],
    output:
        bed = join(outdir, '07_joint_vcf/02_variant_annotations/AM_variants.bed'),
        table = join(outdir, '07_joint_vcf/02_variant_annotations/AlphaMissense_simple_table.tsv')
    shell: """
        grep -v "^#" {input.vcf} | awk 'BEGIN{{OFS="\\t"}} {{print $1, $2-1, $2}}' | sort -k1,1 -k2,2n -u > {output.bed}
        printf "CHROM\\tPOS\\tREF\\tALT\\tgenome\\tuniprot_id\\ttranscript_id\\tprotein_variant\\tam_pathogenicity\\tam_class\\n" > {output.table}
        tabix {input.am_file} -R {output.bed} >> {output.table}
    """

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
        ls x*.vcf | xargs -P {threads} -I {{}} sh -c "$COMMAND -a -c 1 -g {genome_version} -o {params.tmpout}/{{}} {{}}.gz"

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
            --ref-version {funcotator_ref_version} \
            --data-sources-path {params.funcotator_data_path}  \
            --output {output.vcf} \
            --output-file-format VCF
        """

# CADD and AlphaMissense are both optional deleteriousness annotations here (each only
# runs if config asks for it) - see skip_cadd/alphamissense_file comments in setup_germline.smk
def _cadd_input_germline(wildcards):
    return [] if config['skip_cadd'] else rules.subset_columns.output.vcf

def _am_input_germline(wildcards):
    return rules.alphamissense.output.table if config['alphamissense_file'] else []

# combine the annotation from CADD/AlphaMissense and funcotator
rule combine_annoations:
    input:
        vcf_funcotator = rules.funcotator.output.vcf,
        vcf_cadd = _cadd_input_germline,
        vcf_alphamissense = _am_input_germline,
    output:
        vcf = join(outdir, '07_joint_vcf/02_variant_annotations/annotations_combined.vcf'),
    params:
        filter_min_dp = config['filter_min_dp'],
    script: "combine_annotation_germline.R"
