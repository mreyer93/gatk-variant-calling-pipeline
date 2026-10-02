#####################################################################################################
### TUMOR VS NORMAL CALLING #########################################################################
#####################################################################################################
# use empty string for targets if targets file is not specified
if targets != '':
    target_input_list = [targets]
    targets_string = '-L ' + targets
else:
    target_input_list = []
    targets_string = ''

rule HaplotypeCaller:
    input: 
        target_input_list,
        bam = lambda wildcards: get_final_bam(wildcards.sample),
        bai = lambda wildcards: get_final_bai(wildcards.sample),
        ref = REF_FILE,
    output:
        gvcf = join(outdir, '05_haplotypecaller/{sample}.g.vcf.gz')
    threads: 8
    shell: """
        gatk HaplotypeCaller --java-options "-Xmx7g" \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.gvcf} \
            -ERC GVCF \
            {targets_string} \
            --native-pair-hmm-threads {threads} 
    """

rule create_map_file:
    input: 
        expand(join(outdir, '05_haplotypecaller/{sample}.g.vcf.gz'), sample=sample_list)
    output:
        gvcf = join(outdir, '06_GDB/sample_map.txt'),
    threads: 2
    run: 
        with open(output[0], "w") as outf:
            for sample in sample_list:
                gvcf_file = join(outdir, f'05_haplotypecaller/{sample}.g.vcf.gz')
                outf.write(f"{sample}\t{gvcf_file}\n")

rule GenomicsDBImport:
    input: 
        ref = REF_FILE,
        sample_name_map = rules.create_map_file.output
    output:
        GDB = join(outdir, '06_GDB/{chromosome}/FILE'),
    params:
        outdir = join(outdir, '06_GDB/{chromosome}')
    threads: 4
    shell: """
        rm -r {params.outdir}
        gatk --java-options "-Xmx20g" GenomicsDBImport \
            --reference {input.ref} \
            --sample-name-map {input.sample_name_map} \
            --genomicsdb-workspace-path {params.outdir} \
            --reader-threads {threads} \
            --batch-size 32 \
            -L {wildcards.chromosome}
        touch {output}
    """

rule GenotypeGVCFs:
    input: 
        ref = REF_FILE,
        GDB = join(outdir, '06_GDB/{chromosome}/FILE'),
    output:
        vcf = join(outdir, '06_GDB/{chromosome}/{chromosome}.vcf'),
    params:
        workspace = join(outdir, '06_GDB/{chromosome}')
    threads: 4
    shell: """
        gatk GenotypeGVCFs \
            --reference {input.ref} \
            --variant gendb://{params.workspace} \
            --output {output.vcf}
    """

rule mergeVCFs:
    input:
        expand(join(outdir, '06_GDB/{chromosome}/{chromosome}.vcf'), chromosome = chromosome_list)
    output:
        join(outdir, '07_joint_vcf/germline_calls_unfiltered.vcf')
    threads: 1
    params:
        input_vcf_str = ' '.join(['I='+i for i in expand(join(outdir, '06_GDB/{chromosome}/{chromosome}.vcf'), chromosome = chromosome_list)])
    shell: """
        picard MergeVcfs \
            {params.input_vcf_str} \
            O={output}
    """

rule selectSNPs:
    input:
        rules.mergeVCFs.output
    output:
        join(outdir, '07_joint_vcf/germline_calls_SNP.vcf')
    threads:1
    shell: """
        gatk SelectVariants \
            -V {input} \
            -select-type SNP \
            -O {output}
    """

rule hardFilter:
    input:
        rules.selectSNPs.output
    output:
        join(outdir, '07_joint_vcf/germline_calls_SNP_hard_filter.vcf')
    threads:1
    shell: """
        gatk VariantFiltration \
            -V {input} \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "SOR > 3.0" --filter-name "SOR3" \
            -filter "FS > 60.0" --filter-name "FS60" \
            -filter "MQ < 40.0" --filter-name "MQ40" \
            -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O {output}
    """

# Indels (and mixed SNP/indel sites) get GATK's indel hard filters, which differ from the
# SNP ones; without this step every indel was dropped from the filtered calls. Thresholds
# from GATK's "Hard-filtering germline short variants" (article 360035531112).
rule selectIndels:
    input:
        rules.mergeVCFs.output
    output:
        join(outdir, '07_joint_vcf/germline_calls_INDEL.vcf')
    threads:1
    shell: """
        gatk SelectVariants \
            -V {input} \
            -select-type INDEL \
            -select-type MIXED \
            -O {output}
    """

rule hardFilter_indels:
    input:
        rules.selectIndels.output
    output:
        join(outdir, '07_joint_vcf/germline_calls_INDEL_hard_filter.vcf')
    threads:1
    shell: """
        gatk VariantFiltration \
            -V {input} \
            -filter "QD < 2.0" --filter-name "QD2" \
            -filter "QUAL < 30.0" --filter-name "QUAL30" \
            -filter "FS > 200.0" --filter-name "FS200" \
            -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
            -O {output}
    """

# SNPs and indels back together, filter annotations kept
rule mergeFiltered:
    input:
        snps = rules.hardFilter.output,
        indels = rules.hardFilter_indels.output
    output:
        join(outdir, '07_joint_vcf/germline_calls_hard_filter.vcf')
    threads:1
    shell: """
        picard MergeVcfs \
            I={input.snps} \
            I={input.indels} \
            O={output}
    """

rule hardFilter_select:
    input:
        rules.mergeFiltered.output
    output:
        join(outdir, '07_joint_vcf/germline_calls_hard_filter_select.vcf')
    threads:1
    shell: """
        gatk SelectVariants \
            -V {input} \
            --exclude-filtered \
            -O {output}
    """

rule CollectVariantCallingMetrics:
    input:
        vcf = rules.mergeFiltered.output,
        dbsnp = dbsnp_file
    output:
        join(outdir, '08_germline_metrics/done.tmp')
    threads:1
    params:
        outdir = join(outdir, '08_germline_metrics/')
    shell: """
        gatk CollectVariantCallingMetrics \
            -I {input.vcf} \
            --DBSNP  {input.dbsnp} \
            -O {params.outdir}
        touch {output}
    """

rule roh:
    input: 
        vcf = join(outdir, '06_GDB/{chromosome}/{chromosome}.vcf'),
    output:
        roh = join(outdir, '07_roh/{chromosome}/{sample}_roh.txt.gz')
    params:
        roh_intermediate = join(outdir, '07_roh/{chromosome}/{sample}_roh.txt')
    threads: 2
    shell: """
        bcftools roh \
            --samples {wildcards.sample} \
            --threads {threads} \
            -G30 {input.vcf} > \
            {params.roh_intermediate}
        pigz -p {threads} {params.roh_intermediate}
    """

# what's the total amount of ROH that is encompassed in this sample?
# bcftools roh writes one RG line per run of homozygosity (sample, chromosome, start,
# end, length, number of markers, quality) and one ST line per site; the stats come from
# the RG lines, the tool's own segment calls. This used to rebuild segments from the ST
# lines with read_csv(comment='R'), which cut every line at the first capital R (any
# sample name containing one, e.g. NORMAL, crashed it) and dropped any run reaching the
# end of the chromosome.
rule aggregate_roh:
    input:
        roh_files = lambda wildcards: expand(join(outdir, '07_roh/{chromosome}/{sample}_roh.txt.gz'),
                                             chromosome=wildcards.chromosome, sample=sample_list)
    output:
        df = join(outdir, '08_roh_stats/{chromosome}_roh_stats.tsv')
    threads: 1
    run:
        rows = []
        for sample, f in zip(sample_list, input.roh_files):
            with gzip.open(f, 'rt') as fh:
                lengths = [int(line.split('\t')[5]) for line in fh if line.startswith('RG\t')]
            rows.append({'sample': sample, 'roh_number': len(lengths),
                         'roh_total_length': sum(lengths)})
        pd.DataFrame(rows).to_csv(output.df, sep='\t', index=False)