################################################################################
# normalize read names to sample names, which makes all the downstream steps
# much easier!
rule read_symlinks:
    input: 
        fwd = lambda wildcards: read_map[wildcards.sample][0],
        rev = lambda wildcards: read_map[wildcards.sample][1],
    output:
        fwd = join(outdir, '00_read_symlinks/{sample}_R1.fastq.gz'),
        rev = join(outdir, '00_read_symlinks/{sample}_R2.fastq.gz'),
    shell: """
        ln -s $(readlink -f {input[0]}) {output[0]}
        ln -s $(readlink -f {input[1]}) {output[1]}
    """
################################################################################
rule pre_fastqc:
    input:
        fwd = join(outdir, '00_read_symlinks/{sample}_R1.fastq.gz'),
        rev = join(outdir, '00_read_symlinks/{sample}_R2.fastq.gz'),
    output:
        join(outdir,  "00_qc_reports/pre_fastqc/{sample}_R1_fastqc.html"),
        join(outdir,  "00_qc_reports/pre_fastqc/{sample}_R2_fastqc.html")
    params:
        outdir = join(outdir, "00_qc_reports/pre_fastqc/")
    shell: """
        mkdir -p {params.outdir}
        fastqc {input} --outdir {params.outdir}
    """

rule pre_multiqc:
    input: expand(join(outdir, "00_qc_reports/pre_fastqc/{sample}_{read}_fastqc.html"), sample=samples, read=['R1', 'R2'])
    output: report(join(outdir,  "00_qc_reports/pre_multiqc/multiqc_report.html"), category='qc')
    params:
        indir  = join(outdir, "00_qc_reports/pre_fastqc"),
        outdir = join(outdir, "00_qc_reports/pre_multiqc/")
    singularity: "docker://quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
    shell: """
        multiqc --force {params.indir} -o {params.outdir}
    """

################################################################################
rule trim_galore:
    input:
        fwd = join(outdir, '00_read_symlinks/{sample}_R1.fastq.gz'),
        rev = join(outdir, '00_read_symlinks/{sample}_R2.fastq.gz'),
    output:
        fwd = join(outdir, "01_trimmed/{sample}_R1_val_1.fq.gz"), 
        rev = join(outdir, "01_trimmed/{sample}_R2_val_2.fq.gz"),
    threads: 2
    params:
        q_min   = config['trim_galore']['quality'],
        min_len = config['trim_galore']['min_read_length'],
        outdir  = join(outdir, "01_trimmed/")
    log:
        join(outdir, "logs/{sample}_trim.log")
    shell: """
        mkdir -p {params.outdir}
        trim_galore --quality {params.q_min} \
            --length {params.min_len} \
            --cores {threads} \
            --output_dir {params.outdir} \
            --paired {input.fwd} {input.rev}
    """

################################################################################
rule post_fastqc:
    input:
        fwd = join(outdir, "01_trimmed/{sample}_R1_val_1.fq.gz"),
        rev = join(outdir, "01_trimmed/{sample}_R2_val_2.fq.gz")
    output: 
        fwd = join(outdir,  "00_qc_reports/post_fastqc/{sample}_R1_val_1_fastqc.html"),
        rev = join(outdir,  "00_qc_reports/post_fastqc/{sample}_R2_val_2_fastqc.html")
    params:
        outdir = join(outdir, "00_qc_reports/post_fastqc/")
    shell: """
        fastqc {input} --outdir {params.outdir}
    """

rule post_multiqc:
    input: 
        fwd = expand(join(outdir,  "00_qc_reports/post_fastqc/{sample}_R1_val_1_fastqc.html"), sample=samples),
        rev = expand(join(outdir,  "00_qc_reports/post_fastqc/{sample}_R2_val_2_fastqc.html"), sample=samples)
    output: report(join(outdir, "00_qc_reports/post_multiqc/multiqc_report.html"), category='qc')
    params:
        indir  = join(outdir,  "00_qc_reports/post_fastqc"),
        outdir = join(outdir,  "00_qc_reports/post_multiqc/")
    singularity: "docker://quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
    shell: """
        multiqc --force {params.indir} -o {params.outdir}
    """

################################################################################
# after coverage across exons is calculated, run the QC stat script 
# via Rmarkdown
rule coverage_report:
    input:
        join(outdir, "02_align/flagstat_offtarget.txt"),
        expand(join(outdir, "02_align/cov/{sample}.mkdup.bam.w100.cov"), sample=samples)
    output:
        join(outdir, 'coverage_quality_report.pdf')
    params:
        outdir = outdir,
        target_bed = target_bed,
        sample_reads = sample_reads_f
    conda: "envs/rmarkdown.yml"

    script: """
        scripts/primer_check.Rmd
    """
