#####################################################################################################
### CLIENT-FACING SUMMARY REPORT ####################################################################
#####################################################################################################
# Builds a single self-contained report (HTML, plus PDF if enabled) summarising the run:
# QC, variant counts, functional consequences, recurrently mutated genes, a ranked
# shortlist of prioritised variants, and copy-number results. Every figure's underlying
# table is also written to 09_report/summary_tables/ as TSV.
#
# The report adapts to how the pipeline was run - if annotation was skipped
# (config_call_bam_GATK_local.yaml), the gene/consequence/deleteriousness sections are
# omitted and the QC + variant-count sections still render.

report_scripts_dir = workflow.current_basedir


def _report_inputs(wildcards):
    """Dependencies computed at DAG time, matching whatever this run actually produces."""
    deps = []
    if not config['annotation_only']:
        deps += expand(
            join(outdir, '02_variants_reference/01_gatk_variant_calling/filtered/{sample}_filtered.vcf'),
            sample=sample_list)
    if not config.get('skip_annotation', False):
        deps += expand(
            join(outdir, '02_variants_reference/02_variant_annotations/annotations_combined/filtered/{sample}_filtered_annotated.vcf'),
            sample=sample_list)
        deps += expand(
            join(outdir, '02_variants_reference/03_variant_new_scores/final/filtered/{pt}.vcf'),
            pt=unique_pts)
    return deps


rule somatic_report_html:
    input: _report_inputs
    output: join(outdir, '09_report/somatic_report.html')
    log: join(outdir, '09_report/logs/somatic_report_html.log')
    params:
        rmd = join(report_scripts_dir, 'somatic_report.Rmd'),
        scripts_dir = report_scripts_dir,
        format = 'html_document',
        outdir = outdir,
        metadata = config['bam_metadata'],
        project_name = config.get('project_name', 'Somatic variant analysis'),
        genome_version = genome_version,
        filter_min_af = config['filter_min_af'],
        filter_min_dp = config['filter_min_dp'],
        filter_min_cadd = config['filter_min_cadd'],
        filter_gatk = config['filter_GATK'],
        top_n = config.get('report_top_n', 40),
    conda: '../../envs/rmarkdown.yml'
    script: 'render_report.R'


rule somatic_report_pdf:
    input: _report_inputs
    output: join(outdir, '09_report/somatic_report.pdf')
    log: join(outdir, '09_report/logs/somatic_report_pdf.log')
    params:
        rmd = join(report_scripts_dir, 'somatic_report.Rmd'),
        scripts_dir = report_scripts_dir,
        format = 'pdf_document',
        outdir = outdir,
        metadata = config['bam_metadata'],
        project_name = config.get('project_name', 'Somatic variant analysis'),
        genome_version = genome_version,
        filter_min_af = config['filter_min_af'],
        filter_min_dp = config['filter_min_dp'],
        filter_min_cadd = config['filter_min_cadd'],
        filter_gatk = config['filter_GATK'],
        top_n = config.get('report_top_n', 40),
    conda: '../../envs/rmarkdown.yml'
    script: 'render_report.R'
