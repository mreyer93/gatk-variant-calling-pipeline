################################################################################
# Aligns reads to specified reference genome and sorts the output. 
# Only retains and sorts reads where both read mates align to the reference.
rule align_to_ref:
    input:
        bwa_build = rules.build_ref_index.output,
        bwa_index = REF_FILE,
        fwd = rules.trim_galore.output.fwd,
        rev = rules.trim_galore.output.rev
    output:
        withdup=join(outdir, "02_align/bam/{sample}.bam"),
        mkdup=join(outdir, "02_align/bam/{sample}.mkdup.bam"),
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tLB:lib1\tPL:ILLUMINA\tPU:unit1",
    threads: 4
    shell: """
        bwa mem -R '{params.rg}' -t {threads} {input.bwa_index} {input.fwd} {input.rev} \
        | samtools view -b -F12 - \
        | samtools fixmate -m - - \
        | samtools sort -O BAM \
        | tee {output.withdup} \
        | samtools markdup -r - {output.mkdup}
    """

################################################################################
rule count_aligned_reads:
    input:
        expand(join(outdir, "02_align/bam/{sample}.mkdup.bam"), sample=samples)
    output:
        report(join(outdir, '02_align/aligned_counts.txt'), category='alignment')
    shell: """
        for i in {input}; do
            samtools index "$i"
            bn=$(basename "$i" | sed 's/.bam//g')
            echo "$bn\t$(samtools view -c -F 260 $i)" >> {output}
        done
    """

################################################################################
# coverage over exons and over exons +/-100bp
rule coverage_bed_0:
    input:
        bam = rules.align_to_ref.output.mkdup,
        target_bed = target_bed,
        ref = REF_FILE + '.fai',
    output:
        cov = join(outdir, "02_align/cov/{sample}.mkdup.bam.w0.cov"),
    params:
        tmpfile = join(outdir, "02_align/cov/{sample}.mkdup.bam.w0.cov.tmp"),
        read_length = read_length,
        aml_cov = join(outdir, "02_align/out.aml.cov")
    shell: """
        cut -f 1,2 {input.ref} > {output.cov}.genome
        coverageBed -sorted -g {output.cov}.genome -a {input.target_bed} -b {input.bam} > {params.tmpfile}
        awk -v r={params.read_length}  '{{ print $6*r/$8; }}' {params.tmpfile} | paste {params.tmpfile}  - > {output.cov}
        rm {params.tmpfile} {output.cov}.genome
    """

rule coverage_bed_100:
    input:
        bam = rules.align_to_ref.output.mkdup,
        target_bed = target_bed_w100,
        ref = REF_FILE + '.fai',
    output:
        cov = join(outdir, "02_align/cov/{sample}.mkdup.bam.w100.cov"), 
    params:
        tmpfile = join(outdir, "02_align/cov/{sample}.mkdup.bam.w100.cov.tmp"),
        read_length = read_length,
        aml_cov = join(outdir, "02_align/out.aml.cov")
    shell: """
        cut -f 1,2 {input.ref} > {output.cov}.genome
        coverageBed -sorted -g {output.cov}.genome -a {input.target_bed} -b {input.bam} > {params.tmpfile}
        awk -v r={params.read_length}  '{{ print $6*r/$8; }}' {params.tmpfile} | paste {params.tmpfile}  - > {output.cov}
        rm {params.tmpfile} {output.cov}.genome
    """


################################################################################
# flagstats on mkdup.bam file
rule flagstats:
    input:
        expand(join(outdir, "02_align/bam/{sample}.mkdup.bam"), sample=samples)
    output:
        report(join(outdir, "02_align/out_flagstat.txt"), category='alignment')
    shell: """
        for s in {input}; do 
            paste <(echo $s) <(samtools flagstat $s) >> {output}
        done
    """

################################################################################
# identify offtarget and run flagstat on those
rule flagstats_offtarget:
    input:
        bam = join(outdir, "02_align/bam/{sample}.mkdup.bam"),
        target_bed_w1000 = target_bed_w1000,
        ref = REF_FILE + ".fai",
    output:
        bam = join(outdir, "02_align/offtarget/{sample}_off.bam"),
        flagstat = join(outdir, "02_align/offtarget/{sample}_off.bam.flagstat")
    params:
        sample = "{sample}",
        flagstat_out = join(outdir, "02_align/out_flagstat_offtarget.txt")
    shell: """
        cut -f 1,2 {input.ref} > {params.sample}.genome
        intersectBed -sorted -g {params.sample}.genome -v -abam {input.bam} -b {input.target_bed_w1000} > {output.bam}
        echo {output.bam} > {output.flagstat}
        samtools flagstat {output.bam} >> {output.flagstat}
        rm {params.sample}.genome
    """

rule flagstats_aggregate:
    input:
        expand(join(outdir, "02_align/offtarget/{sample}_off.bam.flagstat"), sample=samples)
    output:
        report(join(outdir, "02_align/flagstat_offtarget.txt"), category='alignment')
    shell: """
        cat {input} > {output}
    """

################################################################################
# Rmarkdown QC script
rule primer_check:
    input:
        expand(join(outdir, "02_align/bam/{sample}.mkdup.bam"), sample=samples)
    output: 
        join(outdir, 'primer_check.pdf')
    params: 
        outdir = outdir, 
        sample_reads = config["sample_file"],
        target_bed = config["target_bed"]
    conda:
        '../../envs/rmarkdown.yml'
    script:
        'primer_check.Rmd'
