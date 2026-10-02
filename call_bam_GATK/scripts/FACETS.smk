#####################################################################################################
### FACETS FOR CNA AND LOH ##########################################################################
#####################################################################################################
# FACETS is built from its pinned upstream release by rule facets_build, on every platform,
# inside envs/facets.yml (see that file for why bioconda's snp-pileup and r-facets are not
# used). facets is GPL >= 2, so it is downloaded at build time rather than vendored here.
FACETS_VERSION = '0.6.2'
FACETS_TOOLS = join(outdir, '04_FACETS/tools/facets-' + FACETS_VERSION)

# preProcSample looks up GC content by position in build-specific tables and defaults to
# hg19, so the build has to be passed explicitly for GRCh38 data.
FACETS_GBUILD = {'GRCh38': 'hg38', 'hg38': 'hg38', 'GRCh37': 'hg19', 'hg19': 'hg19'}.get(str(config['genome_version']))

rule facets_build:
    output:
        snp_pileup = join(FACETS_TOOLS, 'bin/snp-pileup'),
        rlib = directory(join(FACETS_TOOLS, 'R'))
    log: join(outdir, '04_FACETS/logs/facets_build.log')
    params:
        url = 'https://github.com/mskcc/facets/archive/refs/tags/v' + FACETS_VERSION + '.tar.gz'
    conda: "../../envs/facets.yml"
    shell: """
        set -euo pipefail
        src=$(mktemp -d)
        curl -sSfL {params.url} | tar -xz -C "$src" --strip-components 1
        mkdir -p "$(dirname {output.snp_pileup})" {output.rlib}
        # argp is part of glibc on Linux; macOS needs argp-standalone's libargp
        ARGP=""
        if [ "$(uname -s)" = "Darwin" ]; then ARGP="-largp"; fi
        "$CXX" $CXXFLAGS $CPPFLAGS -std=c++11 -I"$CONDA_PREFIX/include" \
            "$src/inst/extcode/snp-pileup.cpp" $LDFLAGS -L"$CONDA_PREFIX/lib" -lhts $ARGP \
            -Wl,-rpath,"$CONDA_PREFIX/lib" -o {output.snp_pileup} > {log} 2>&1
        R CMD INSTALL --library={output.rlib} "$src" >> {log} 2>&1
        rm -rf "$src"
    """

# takes in normal and tumor samples for the same patient
rule facets_snp_pileup:
    input:
        bams = lambda wildcards: [get_final_bam(sample) for sample in pt_tp_to_samples[wildcards.patient_tp]],
        dbsnp_common_file = dbsnp_common_file,
        snp_pileup = rules.facets_build.output.snp_pileup
    output:
        csv = join(outdir, '04_FACETS/{patient_tp}.csv.gz')
    params:
        bam_string = lambda wildcards: ' '.join(pt_tp_to_final_bams[wildcards.patient_tp]),
    conda: "../../envs/facets.yml"
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
        {input.snp_pileup} -g -q15 -Q20 -P100 -r25,0 -v -d 10000 \
            {input.dbsnp_common_file} {output} {params.bam_string}
    """

rule facets_snp_plot:
    input:
        csv = rules.facets_snp_pileup.output.csv,
        rlib = rules.facets_build.output.rlib
    params:
        patient_tp = lambda wildcards: wildcards.patient_tp,
        normal_samples = lambda wildcards: (n_samples_pt_tp[wildcards.patient_tp]),
        tumor_samples = lambda wildcards: (t_samples_pt_tp[wildcards.patient_tp]),
        gbuild = FACETS_GBUILD,
        # one SNP per window of roughly the insert size: 250 for exomes (the FACETS default),
        # ~150 for targeted panels or short inserts, ~500 for WGS (FACETS issues #61, #81)
        snp_nbhd = config.get('facets_snp_nbhd', 250)
    output:
        pdf = join(outdir, '04_FACETS/{patient_tp}.pdf'),
        txt = join(outdir, '04_FACETS/{patient_tp}_purity.txt'),
        cncf = join(outdir, '04_FACETS/{patient_tp}_cncf.tsv')
    conda: "../../envs/facets.yml"
    script: "facets_plotting.R"
