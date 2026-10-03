#!/usr/bin/env bash
# Checks that scattering the callers over interval shards does not change any call.
#
# Runs the somatic and germline pipelines twice each on the smoke-test data
# (test/run_test.sh prepares it), over the same tiled target region: once as a single
# shard and once split into four. Every call (position, alleles, FILTER, INFO and
# genotypes) must be identical between the two: in the three somatic call sets (each
# sample against the reference, and tumour vs normal) and in the germline joint calls.
#
# Exact equality holds here because no position in this data has more than 50 reads
# starting at it. Where one does, GATK keeps 50 chosen with a random sequence shared across
# the run, so that site's result depends on what was processed before, scattered or not; on
# a real exome this moved one PASS call per call set (see manual/somatic.md).
#
# Usage: ./test/test_scatter.sh

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

CORES="${CORES:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
D=test/data
REF="$D/reference"

if [[ ! -s "$D/config_test.yaml" || ! -s "$D/bam/TUMOUR.bam" ]]; then
    echo "==> Preparing the shared test data (test/run_test.sh --dry-run)"
    ./test/run_test.sh --dry-run > /dev/null
fi

# Ten 1 kb tiles with 100 bp gaps across the region the reads cover (chr1:131,000-142,000),
# so the shards have calls on both sides of each boundary.
awk 'BEGIN{OFS="\t"; for (i = 0; i < 10; i++) {s = 131000 + i * 1100; print "chr1", s, s + 1000, "+", "tile_" i}}' \
    > "$REF/tiles.bed"

fail=0
[[ -s "$D/config_test_germline.yaml" ]] || ./test/run_test_germline.sh --dry-run > /dev/null

for n in 1 4; do
    out="$D/results_scatter$n"
    sed -e "s|^output_directory:.*|output_directory: \"$out\"|" \
        -e "s|^target_bed_w100:.*|target_bed_w100: \"$REF/tiles.bed\"|" \
        -e "s|^make_report:.*|make_report: False|" \
        "$D/config_test.yaml" > "$D/config_scatter$n.yaml"
    sed -e "s|^output_directory:.*|output_directory: \"$out/germline\"|" \
        -e "s|^targets:.*|targets: \"$REF/tiles.bed\"|" \
        "$D/config_test_germline.yaml" > "$D/config_scatter${n}_germline.yaml"
    echo "==> scatter_count: $n"
    snakemake -s call_bam_GATK/call_bam_GATK.snakefile --configfile "$D/config_scatter$n.yaml" \
        --config scatter_count=$n --cores "$CORES" -q > "$D/scatter$n.log" 2>&1 \
        || { echo "  somatic pipeline failed, see $D/scatter$n.log" >&2; exit 1; }
    snakemake -s call_bam_GATK/call_bam_GATK_germline.snakefile --configfile "$D/config_scatter${n}_germline.yaml" \
        --config scatter_count=$n --cores "$CORES" -q > "$D/scatter${n}_germline.log" 2>&1 \
        || { echo "  germline pipeline failed, see $D/scatter${n}_germline.log" >&2; exit 1; }
    echo "  $(ls "$out/00_scatter" | wc -l | tr -d ' ') shard(s)"
done

for rel in 02_variants_reference/01_gatk_variant_calling/filtered/TUMOUR_filtered.vcf \
           02_variants_reference/01_gatk_variant_calling/filtered/NORMAL_filtered.vcf \
           03_variants_TvN/01_gatk_variant_calling/filtered/PATIENT_A_Dx_filtered.vcf \
           germline/07_joint_vcf/germline_calls_hard_filter.vcf; do
    a=$(grep -v '^#' "$D/results_scatter1/$rel" | sort)
    b=$(grep -v '^#' "$D/results_scatter4/$rel" | sort)
    n=$(printf '%s\n' "$a" | grep -c . || true)
    if [[ "$a" == "$b" ]]; then
        echo "  OK   $(basename "$rel"): $n calls, identical"
    else
        echo "  DIFF $(basename "$rel")"; diff <(printf '%s\n' "$a") <(printf '%s\n' "$b") | head -10; fail=1
    fi
done
echo
if [[ $fail -eq 0 ]]; then echo "Scatter test PASSED."; else echo "Scatter test FAILED." >&2; exit 1; fi
