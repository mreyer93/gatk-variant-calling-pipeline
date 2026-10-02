#!/usr/bin/env bash
# End-to-end smoke test for the read-processing pipeline (process_reads): FastQC ->
# Trim Galore -> bwa mem -> duplicate removal -> coverage and off-target QC -> the
# coverage report (PDF).
#
# Uses the same public data as test/run_test.sh (the nf-core/sarek tiny pair), with each
# sample's lanes concatenated into one FASTQ pair, since process_reads takes one pair per
# sample. The capture panel is a made-up one sized to the data: three 1 kb targets inside
# the region the reads cover, and one target with no reads, so the report's low-coverage
# and off-target sections have something to show.
#
# Usage:
#   ./test/run_test_processing.sh              # prepare the shared test data if needed, then run
#   ./test/run_test_processing.sh --dry-run    # build the DAG only, run nothing

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

DRY_RUN=""
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN="-n"

# The coverage report renders in a Snakemake-managed conda env (envs/rmarkdown.yml).
CONDA_FLAG="--use-conda"
command -v mamba >/dev/null 2>&1 || CONDA_FLAG="$CONDA_FLAG --conda-frontend conda"
[[ "${USE_CONDA:-1}" == "0" ]] && CONDA_FLAG=""
# R cannot run from a path containing a space; keep the envs somewhere that has none.
if [[ -n "$CONDA_FLAG" && "$PWD" == *" "* ]]; then
    CONDA_FLAG="$CONDA_FLAG --conda-prefix ${SNAKEMAKE_CONDA_PREFIX:-$HOME/.cache/snakemake-conda}"
fi

# Use every available core by default; set CORES to cap it.
CORES="${CORES:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"

D=test/data
REF="$D/reference"; FQ="$D/fastq"; OUT="$D/results_processing"

# The reference (bwa-indexed) and FASTQs are prepared by the somatic smoke test, whose
# download and alignment steps run even in --dry-run mode.
if [[ ! -s "$REF/genome.chr.fasta.bwt" || ! -s "$REF/genome.chr.fasta.fai" ]] \
        || ! ls "$FQ"/tiny_n_L*_R1_xxx.fastq.gz >/dev/null 2>&1; then
    echo "==> Preparing the shared test data (test/run_test.sh --dry-run)"
    ./test/run_test.sh --dry-run > /dev/null
fi

echo "==> Concatenating each sample's lanes"
mkdir -p "$D/fastq_merged"
for s in n:NORMAL t:TUMOUR; do
    for r in R1 R2; do
        # gzip files concatenate into a valid gzip file
        cat "$FQ"/tiny_${s%%:*}_L*_${r}_xxx.fastq.gz > "$D/fastq_merged/${s##*:}_${r}.fastq.gz"
    done
done
printf 'NORMAL\t%s/NORMAL_R1.fastq.gz\t%s/NORMAL_R2.fastq.gz\nTUMOUR\t%s/TUMOUR_R1.fastq.gz\t%s/TUMOUR_R2.fastq.gz\n' \
    "$D/fastq_merged" "$D/fastq_merged" "$D/fastq_merged" "$D/fastq_merged" > "$D/sample_reads_test.tsv"

echo "==> Writing the test panel and config"
# chr, start, end, strand, name; the reads fall in chr1:131,000-142,000
printf 'chr1\t50000\t51000\t+\ttarget_1\nchr1\t132500\t133500\t+\ttarget_2\nchr1\t136000\t137000\t+\ttarget_3\nchr1\t139000\t140000\t+\ttarget_4\n' \
    > "$REF/panel_w0.bed"
for w in 100 1000; do
    awk -v w=$w 'BEGIN{OFS="\t"} {s=$2-w; if (s<0) s=0; print $1, s, $3+w, $4, $5}' "$REF/panel_w0.bed" \
        > "$REF/panel_w$w.bed"
done

cat > "$D/config_test_processing.yaml" <<EOF
output_directory: "$OUT"
sample_file: "$D/sample_reads_test.tsv"
read_length: 101
trim_galore:
  quality: 20
  min_read_length: 30
reference_file: "$REF/genome.chr.fasta"
target_bed: "$REF/panel_w0.bed"
target_bed_w100: "$REF/panel_w100.bed"
target_bed_w1000: "$REF/panel_w1000.bed"
EOF

echo "==> Running pipeline"
snakemake -s process_reads/processing.snakefile \
    --configfile "$D/config_test_processing.yaml" --cores "$CORES" $CONDA_FLAG $DRY_RUN

if [[ -z "$DRY_RUN" ]]; then
    echo; echo "==> Checking expected outputs"
    fail=0
    for f in "$OUT/00_qc_reports/pre_multiqc/multiqc_report.html" \
             "$OUT/00_qc_reports/post_multiqc/multiqc_report.html" \
             "$OUT/01_trimmed/NORMAL_R1_val_1.fq.gz" \
             "$OUT/02_align/bam/NORMAL.mkdup.bam" "$OUT/02_align/bam/TUMOUR.mkdup.bam" \
             "$OUT/02_align/aligned_counts.txt" "$OUT/02_align/flagstat_offtarget.txt" \
             "$OUT/02_align/cov/NORMAL.mkdup.bam.w100.cov" "$OUT/primer_check.pdf"; do
        if [[ -s "$f" ]]; then echo "  OK   $f"; else echo "  MISS $f"; fail=1; fi
    done
    # Counts: a tab-separated row per sample; duplicate removal only ever removes reads;
    # some reads fall off-target (outside the targets +/- 1 kb) but most do not.
    if [[ $fail -eq 0 ]]; then
        if python3 - "$OUT" <<'PY'
import subprocess, sys
out = sys.argv[1]
problems = []
counts = [l.split("\t") for l in open(f"{out}/02_align/aligned_counts.txt").read().splitlines()]
if sorted(c[0] for c in counts) != ["NORMAL.mkdup", "TUMOUR.mkdup"] or any(len(c) != 2 for c in counts):
    problems.append(f"aligned_counts.txt is not one tab-separated row per sample: {counts}")
def n(bam):
    return int(subprocess.check_output(["samtools", "view", "-c", "-F", "0x904", bam]))
for s in ["NORMAL", "TUMOUR"]:
    before, after = n(f"{out}/02_align/bam/{s}.bam"), n(f"{out}/02_align/bam/{s}.mkdup.bam")
    off = n(f"{out}/02_align/offtarget/{s}_off.bam")
    if not 0 < after <= before:
        problems.append(f"{s}: {before} reads before duplicate removal, {after} after")
    elif not 0 < off < after / 2:
        problems.append(f"{s}: {off} of {after} reads off-target")
    else:
        print(f"  OK   {s}: {before} reads aligned, {after} after duplicate removal, "
              f"{off} ({100 * off / after:.1f}%) off-target")
for p in problems:
    print("  CHECK " + p)
sys.exit(1 if problems else 0)
PY
        then :; else fail=1; fi
    fi
    echo
    if [[ $fail -eq 0 ]]; then
        echo "Processing smoke test PASSED."
    else
        echo "Processing smoke test FAILED - see the Snakemake log above." >&2; exit 1
    fi
fi
