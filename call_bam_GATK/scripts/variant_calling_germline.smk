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
        bai = lambda wildcards: get_final_bam(wildcards.sample)[:-1] + 'i',
        # bai = rules.bam_index.output,
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

rule hardFilter_select:
    input:
        rules.hardFilter.output
    output:
        join(outdir, '07_joint_vcf/germline_calls_SNP_hard_filter_select.vcf')
    threads:1
    shell: """
        gatk SelectVariants \
            -V {input} \
            --exclude-filtered \
            -O {output}
    """

rule CollectVariantCallingMetrics:
    input:
        vcf = rules.hardFilter.output,
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
def calculate_roh_df(df):
    current_state=0
    counter=0
    roh_starts = []
    roh_ends = []
    for roh, pos in zip(df.roh, df.pos):
        if current_state==0:
            if roh ==1:
                counter +=1 
                roh_starts.append(pos)
                current_state=1
        else: 
            if roh==0:
                current_state=0
                roh_ends.append(pos)
    roh_calcs = pd.DataFrame({'start': roh_starts[0:len(roh_ends)], 
                          'end' : roh_ends[0:len(roh_ends)]})
    roh_calcs['len'] = [a-b for a,b in zip(roh_calcs['end'], roh_calcs['start'])]    
    return(roh_calcs)

rule aggregate_roh:
    input: 
        roh_files = expand(join(outdir, '07_roh/{chromosome}/{sample}_roh.txt.gz'),chromosome=chromosome_list, sample=sample_list)
    output:
        df = join(outdir, '08_roh_stats/{chromosome}_roh_stats.tsv')
    params:
        basedir = join(outdir, '07_roh/{chromosome}')
    threads: 1
    run: 
        basedir = params.basedir
        roh_calcs_dict = {}
        total_roh_length_dict = {}
        roh_number_dict = {}
        for sample in sample_list:
            f = join(basedir, f"{sample}_roh.txt.gz")
            df = pd.read_csv(f, sep='\t', skiprows=5, comment='R', header=None)
            df.columns = ['class', 'sample', 'chr', 'pos', 'roh', 'prob']
            roh_calcs = calculate_roh_df(df)
            total_roh_length = np.sum(roh_calcs['len'])
            roh_number = roh_calcs.shape[0]    
            total_roh_length_dict[sample] = total_roh_length
            roh_number_dict[sample] = roh_number
        roh_stats_df = pd.DataFrame({'sample': sample_list,
                                'roh_number': [roh_number_dict[s] for s in sample_list],
                                'roh_total_length': [total_roh_length_dict[s] for s in sample_list]})
        roh_stats_df.index = roh_stats_df['sample']
        roh_stats_df.to_csv(output[0], sep='\t')