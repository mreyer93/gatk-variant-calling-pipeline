#####################################################################################################
### FACETS FOR CNA AND LOH ##########################################################################
#####################################################################################################
# snp-pileup comes from the `snp-pileup` conda package (envs/environment.yml) rather
# than a vendored binary - install it into the pipeline's conda env so it's on PATH.

# takes in normal and tumor samples for the same patient
rule facets_snp_pileup:
    input:
        bams = lambda wildcards: [get_final_bam(sample) for sample in pt_tp_to_samples[wildcards.patient_tp]],
        dbsnp_common_file = dbsnp_common_file
    output:
        csv = join(outdir, '04_FACETS/{patient_tp}.csv.gz')
    params:
        bam_string = lambda wildcards: ' '.join(pt_tp_to_final_bams[wildcards.patient_tp]),
    shell: """
        set +u
        # options used
        # -g gzip
        # -q min mapping quality
        # -Q min base quality
        # -P pseudo-snps, inert a blank record every interval if no SNPs
        # -r min read count for position to be output in normal, tumor
        # -v verbose output
        # -d max depth
        echo {params.bam_string}
        snp-pileup -g -q15 -Q20 -P100 -r25,0 -v -d 10000 \
            {input.dbsnp_common_file} {output} {params.bam_string}
    """

rule facets_snp_plot:
    input:
        csv = rules.facets_snp_pileup.output
    params:
        patient_tp = lambda wildcards: wildcards.patient_tp,
        normal_samples = lambda wildcards: (n_samples_pt_tp[wildcards.patient_tp]),
        tumor_samples = lambda wildcards: (t_samples_pt_tp[wildcards.patient_tp])
    output:
        pdf = join(outdir, '04_FACETS/{patient_tp}.pdf'),
        txt = join(outdir, '04_FACETS/{patient_tp}_purity.txt')
    script: "facets_plotting.R"
