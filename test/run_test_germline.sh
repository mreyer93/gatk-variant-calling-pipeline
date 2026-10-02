#!/usr/bin/env bash
# End-to-end smoke test for the germline pipeline, on the same public data as
# test/run_test.sh: the nf-core/sarek tiny pair, aligned with bwa mem, treated here as two
# germline samples. Runs read groups -> BQSR -> HaplotypeCaller (GVCF) ->
# GenomicsDBImport -> joint genotyping -> SNP and indel hard filters -> runs of
# homozygosity.
#
# Annotation (Funcotator, CADD, AlphaMissense) is switched off: those data sources are
# gigabytes and are not needed to demonstrate calling. See manual/germline.md.
#
# This checks that the pipeline runs and that its outputs are well formed. It does not
# measure accuracy: the tiny pair has no truth set, and several sites look simulated
# (a NORMAL heterozygote at ~30% alt with zero alt reads in TUMOUR at similar depth), so
# genotype agreement between the two samples is not a meaningful target here.
#
# Usage:
#   ./test/run_test_germline.sh              # prepare the shared test data if needed, then run
#   ./test/run_test_germline.sh --dry-run    # build the DAG only, run nothing

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

DRY_RUN=""
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN="-n"

# Use every available core by default; set CORES to cap it.
CORES="${CORES:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

D=test/data
REF="$D/reference"; BAM="$D/bam"; OUT="$D/results_germline"

# The reference, resources and aligned BAMs are prepared by the somatic smoke test, whose
# download and alignment steps run even in --dry-run mode.
if [[ ! -s "$BAM/TUMOUR.bam" || ! -s "$BAM/NORMAL.bam" || ! -s "$REF/targets.bed" ]]; then
    echo "==> Preparing the shared test data (test/run_test.sh --dry-run)"
    ./test/run_test.sh --dry-run > /dev/null
fi

echo "==> Writing sample sheet and config"
printf 'sample\tbamfile\nNORMAL\t%s/NORMAL.bam\nTUMOUR\t%s/TUMOUR.bam\n' "$BAM" "$BAM" \
    > "$D/bam_metadata_germline.tsv"

cat > "$D/config_test_germline.yaml" <<EOF
project_name: "Germline calling demo (nf-core/sarek tiny pair)"
output_directory: "$OUT"
annotation_only: False
bam_metadata: "$D/bam_metadata_germline.tsv"
filter_min_dp: 5
skip_annotation: True
chromosomes: ["chr1"]   # the tiny reference's reads all map to chr1
targets: "$REF/targets.bed"
genome_version: "GRCh37"
REF_FILE: "$REF/genome.chr.fasta"
funcotator_data_path: "not_used"
dbsnp_file: "$REF/dbsnp_138.chr.vcf.gz"
dbsnp_common_file: "$REF/dbsnp_138.chr.vcf.gz"
gnomad_file: "$REF/gnomAD.chr.vcf.gz"
EOF

echo "==> Running pipeline"
snakemake -s call_bam_GATK/call_bam_GATK_germline.snakefile \
    --configfile "$D/config_test_germline.yaml" --cores "$CORES" $DRY_RUN

if [[ -z "$DRY_RUN" ]]; then
    echo; echo "==> Checking expected outputs"
    fail=0
    J="$OUT/07_joint_vcf"
    for f in "$OUT/05_haplotypecaller/NORMAL.g.vcf.gz" "$OUT/05_haplotypecaller/TUMOUR.g.vcf.gz" \
             "$OUT/06_GDB/chr1/chr1.vcf" \
             "$J/germline_calls_unfiltered.vcf" "$J/germline_calls_hard_filter.vcf" \
             "$J/germline_calls_hard_filter_select.vcf" \
             "$OUT/07_roh/chr1/NORMAL_roh.txt.gz" "$OUT/08_roh_stats/chr1_roh_stats.tsv"; do
        if [[ -s "$f" ]]; then echo "  OK   $f"; else echo "  MISS $f"; fail=1; fi
    done
    # Structure: both samples genotyped, every unfiltered call accounted for by the SNP
    # and indel branches, only PASS records in the selected set, both samples in the ROH
    # table. Then the numbers, for information.
    if [[ $fail -eq 0 ]]; then
        if python3 - "$J" "$OUT/08_roh_stats/chr1_roh_stats.tsv" <<'PY'
import sys
j, roh = sys.argv[1], sys.argv[2]
def records(name):
    rows, samples = [], []
    for line in open(f"{j}/{name}.vcf"):
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
        elif not line.startswith("#"):
            rows.append(line.rstrip("\n").split("\t"))
    return rows, samples
raw, samples = records("germline_calls_unfiltered")
snp, _ = records("germline_calls_SNP")
indel, _ = records("germline_calls_INDEL")
filt, _ = records("germline_calls_hard_filter")
sel, _ = records("germline_calls_hard_filter_select")
problems = []
if sorted(samples) != ["NORMAL", "TUMOUR"]:
    problems.append(f"samples in the joint VCF are {samples}")
if not raw:
    problems.append("no variants called")
if len(snp) + len(indel) != len(raw) or len(filt) != len(raw):
    problems.append(f"{len(raw)} raw calls, but {len(snp)} SNPs + {len(indel)} indels, "
                    f"{len(filt)} after filtering")
if any(r[6] != "PASS" for r in sel):
    problems.append("non-PASS records in the selected calls")
roh_samples = sorted(l.split("\t")[0] for l in open(roh).read().splitlines()[1:])
if roh_samples != ["NORMAL", "TUMOUR"]:
    problems.append(f"ROH table lists {roh_samples}")
for p in problems:
    print("  CHECK " + p)
if not problems:
    print(f"  OK   {len(raw)} joint calls ({len(snp)} SNPs, {len(indel)} indels); "
          f"{len(sel)} pass the hard filters")
sys.exit(1 if problems else 0)
PY
        then :; else fail=1; fi
    fi
    echo
    if [[ $fail -eq 0 ]]; then
        echo "Germline smoke test PASSED."
    else
        echo "Germline smoke test FAILED - see the Snakemake log above." >&2; exit 1
    fi
fi
