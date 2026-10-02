#!/usr/bin/env bash
# FACETS check. The somatic smoke test cannot exercise FACETS: its tumour/normal pair
# covers a few kilobases, far too few heterozygous SNPs to fit copy number. This checks
# the two halves separately instead:
#
#   1. snp-pileup, as built by rule facets_build, reports the same reference and
#      alternate read counts as samtools mpileup with matching filters, at every dbSNP
#      site the smoke-test BAMs cover.
#   2. facets_plotting.R reproduces the worked example that ships with FACETS (a TCGA
#      stomach exome; the vignette reports purity 0.892 and ploidy 2.07) and writes the
#      segment table with its LOH calls.
#
# Run ./test/run_test.sh first: this reuses its BAMs and reference. Needs the
# gatk-pipeline environment active (snakemake, samtools, python3).
#
# Usage:
#   ./test/test_facets.sh

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

D=test/data
CFG=$D/config_test.yaml
TOOLS=$D/results/04_FACETS/tools/facets-0.6.2
[[ -s "$D/bam/NORMAL.bam" && -s "$D/bam/TUMOUR.bam" ]] || { echo "Run ./test/run_test.sh first." >&2; exit 1; }

SMK=(snakemake -s call_bam_GATK/call_bam_GATK.snakefile --configfile "$CFG" --use-conda --cores 4)
command -v mamba >/dev/null 2>&1 || SMK+=(--conda-frontend conda)
# R cannot run from a path containing a space (see run_test.sh)
[[ "$PWD" == *" "* ]] && SMK+=(--conda-prefix "${SNAKEMAKE_CONDA_PREFIX:-$HOME/.cache/snakemake-conda}")

echo "==> Building FACETS (rule facets_build)"
"${SMK[@]}" "$TOOLS/bin/snp-pileup" > /dev/null 2>&1
ENV=$("${SMK[@]}" --list-conda-envs "$TOOLS/bin/snp-pileup" 2>/dev/null | awk -F'\t' '$1 == "envs/facets.yml" {print $3}')
[[ -d "$ENV" ]] || { echo "FAILED: could not find the envs/facets.yml conda env" >&2; exit 1; }
echo

W=$(mktemp -d)
trap 'rm -rf "$W"' EXIT

echo "==> snp-pileup against samtools mpileup"
# snp-pileup drops duplicate, secondary, QC-failed, MAPQ 0 and improperly paired reads and
# counts an overlapping read pair once; samtools mpileup does the same by default. -B
# because snp-pileup applies no base alignment quality.
"$TOOLS/bin/snp-pileup" -q15 -Q20 -r1,0 -d 10000 "$D/reference/dbsnp_138.chr.vcf.gz" "$W/sp.csv" \
    "$D/bam/NORMAL.bam" "$D/bam/TUMOUR.bam" > /dev/null
tail -n +2 "$W/sp.csv" | awk -F, '{print $1"\t"$2}' > "$W/sites.txt"
samtools mpileup -B -q 15 -Q 20 -d 10000 -l "$W/sites.txt" \
    "$D/bam/NORMAL.bam" "$D/bam/TUMOUR.bam" 2>/dev/null > "$W/mpileup.txt"
python3 - "$W/sp.csv" "$W/mpileup.txt" <<'EOF'
import csv, re, sys

def ref_alt(bases, ref, alt):
    s = re.sub(r'\^.', '', bases).replace('$', '')   # read-start (with MAPQ) and read-end marks
    kept, i = [], 0
    while i < len(s):
        if s[i] in '+-':                              # indel: +3ACG / -2TT, skip the sequence
            m = re.match(r'[+-](\d+)', s[i:])
            i += len(m.group(0)) + int(m.group(1))
            continue
        kept.append(s[i].upper())
        i += 1
    return kept.count(ref), kept.count(alt)

sp = {(r['Chromosome'], int(r['Position'])): r for r in csv.DictReader(open(sys.argv[1]))}
n = bad = 0
for line in open(sys.argv[2]):
    x = line.rstrip('\n').split('\t')
    r = sp.get((x[0], int(x[1])))
    if r is None:
        continue
    n += 1
    for f, col in ((1, 4), (2, 7)):
        if ref_alt(x[col], r['Ref'], r['Alt']) != (int(r[f'File{f}R']), int(r[f'File{f}A'])):
            bad += 1
            print(f"  MISMATCH {x[0]}:{x[1]} file {f}", file=sys.stderr)
if n == 0 or n != len(sp) or bad:
    sys.exit(f"FAILED: {n} of {len(sp)} sites compared, {bad} ref/alt mismatches")
print(f"  PASS  ref/alt counts identical at all {n} sites, normal and tumour")
EOF
echo

echo "==> FACETS fit on the stomach example"
cat > "$W/fit.R" <<'EOF'
args <- commandArgs(trailingOnly = TRUE)
setClass("SnakemakeMock", representation(input = "list", output = "list", params = "list"))
snakemake <- new("SnakemakeMock",
    input  = list(csv = file.path(args[2], "facets", "extdata", "stomach.csv.gz"), rlib = args[2]),
    output = list(pdf = file.path(args[3], "stomach.pdf"), txt = file.path(args[3], "purity.txt"),
                  cncf = file.path(args[3], "cncf.tsv")),
    params = list(patient_tp = "stomach", normal_samples = "N", tumor_samples = "T",
                  gbuild = "hg19"))     # the example is a TCGA hg19 exome
invisible(capture.output(source(args[1])))
pp <- read.delim(snakemake@output[["txt"]])
cncf <- read.delim(snakemake@output[["cncf"]])
cat(sprintf("  purity %.4f, ploidy %.4f, %d segments, %d LOH (%d copy-neutral)\n",
            pp$purity, pp$ploidy, nrow(cncf), sum(cncf$loh), sum(cncf$cn_neutral_loh)))
ok <- abs(pp$purity - 0.892) < 0.002 && abs(pp$ploidy - 2.07) < 0.01 &&
      nrow(cncf) > 0 && sum(cncf$loh) > 0 && file.exists(snakemake@output[["pdf"]])
if (!ok) stop("FAILED: does not match the FACETS vignette (purity 0.892, ploidy 2.07)")
cat("  PASS  matches the vignette, and the segment table carries LOH calls\n")
EOF
"$ENV/bin/Rscript" --vanilla "$W/fit.R" call_bam_GATK/scripts/facets_plotting.R "$TOOLS/R" "$W" \
    2> >(grep -v -E 'xcrun|built under' >&2)
echo
echo "FACETS test PASSED"
