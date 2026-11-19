#!/usr/bin/env bash
# Backtests FilterMutectCalls settings on a finished SEQC2 HCC1395 run (see manual/benchmark.md).
#
# Re-filters the run's unfiltered tumour-vs-normal calls three ways and scores each with
# evaluate_seqc2.R:
#   as_run            the pipeline's own filtered VCF (contamination table, orientation priors)
#   segmentation      + --tumor-segmentation from CalculateContamination (GATK best practice)
#   no_contamination  without the contamination table
# FilterMutectCalls takes about 10 s per variant on an exome; the scoring about 20 s.
#
# Usage, from the repository root in the pipeline environment:
#   benchmark/seqc2/backtest_filters.sh <run output_directory> <reference.fasta> <seqc2 dir> <results dir>

set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <run output_directory> <reference.fasta> <seqc2 dir> <results dir>" >&2
    exit 1
fi
run=$1; ref=$2; seqc2=$3; out=$4
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$out/segmentation" "$out/no_contamination"

u="$run/03_variants_TvN/01_gatk_variant_calling/unfiltered"
pileups="$run/01_prepare_bam/GetPileupSummaries"
contam="$run/01_prepare_bam/CalculateContamination/TUMOUR_contam.table"
common=(-R "$ref" -V "$u/HCC1395_Dx_unfiltered.vcf" -stats "$u/HCC1395_Dx_unfiltered.vcf.stats"
        --ob-priors "$u/HCC1395_Dx_readOrientation_model.tar.gz")

echo "==> CalculateContamination with tumour segmentation"
gatk CalculateContamination -I "$pileups/TUMOUR_GetPileupSummaries.table" \
    -matched "$pileups/NORMAL_GetPileupSummaries.table" \
    -O "$out/TUMOUR_contam.table" --tumor-segmentation "$out/TUMOUR_segments.table" > "$out/filter.log" 2>&1
gatk CalculateContamination -I "$pileups/NORMAL_GetPileupSummaries.table" \
    -O "$out/NORMAL_contam.table" >> "$out/filter.log" 2>&1
echo "    tumour contamination: $(awk 'NR==2 {printf "%.4f +/- %.4f", $2, $3}' "$out/TUMOUR_contam.table")"
echo "    normal contamination: $(awk 'NR==2 {printf "%.4f +/- %.4f", $2, $3}' "$out/NORMAL_contam.table")"

echo "==> FilterMutectCalls, two variants"
gatk FilterMutectCalls "${common[@]}" --contamination-table "$contam" \
    --tumor-segmentation "$out/TUMOUR_segments.table" \
    -O "$out/segmentation/HCC1395_Dx_filtered.vcf" >> "$out/filter.log" 2>&1
gatk FilterMutectCalls "${common[@]}" \
    -O "$out/no_contamination/HCC1395_Dx_filtered.vcf" >> "$out/filter.log" 2>&1

echo "==> Scoring"
Rscript "$here/evaluate_seqc2.R" "$run" "$ref" "$seqc2" "$out/eval_as_run" > "$out/eval_as_run.log" 2>&1
for v in segmentation no_contamination; do
    Rscript "$here/evaluate_seqc2.R" "$run" "$ref" "$seqc2" "$out/eval_$v" "$out/$v/HCC1395_Dx_filtered.vcf" \
        > "$out/eval_$v.log" 2>&1
done

Rscript -e '
out <- commandArgs(trailingOnly = TRUE)[1]
d <- do.call(rbind, lapply(c("as_run", "segmentation", "no_contamination"), function(v) {
    s <- read.delim(file.path(out, paste0("eval_", v), "scores.tsv"))
    s <- s[s$callset == "pass", ]; s$filtering <- v; s }))
d <- d[order(d$region, d$type), c("filtering", "region", "type", "calls", "TP", "FP", "FN", "precision", "recall", "F1")]
d[, 8:10] <- round(d[, 8:10], 3)
write.table(d, file.path(out, "backtest_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(d, row.names = FALSE)
# False positives that appear only without the contamination table: a contaminating
# genome would show up as its germline variants, mostly common in gnomAD (POPAF < 2)
fp <- function(v) { f <- read.delim(file.path(out, paste0("eval_", v), "false_positives.tsv"))
    f[f$callset == "pass", ] }
a <- fp("as_run"); b <- fp("no_contamination")
new <- b[!paste(b$CHROM, b$POS, b$REF, b$ALT) %in% paste(a$CHROM, a$POS, a$REF, a$ALT), ]
cat(sprintf("\nFalse positives only without the contamination table: %d (%d SNVs); median tumour AF %.3f\n",
            nrow(new), sum(new$type == "SNV"), median(new$T_AF)))
cat(sprintf("  gnomAD AF >= 1%%: %d; rarer: %d; absent from gnomAD: %d\n",
            sum(new$POPAF < 2), sum(new$POPAF >= 2 & new$POPAF < 6), sum(new$POPAF >= 6)))' "$out"

####
# ==> CalculateContamination with tumour segmentation
#     tumour contamination: 0.0402 +/- 0.0011
#     normal contamination: 0.0003 +/- 0.0001
# ==> FilterMutectCalls, two variants
# ==> Scoring
#         filtering    region  type calls   TP FP  FN precision recall    F1
#            as_run      eval   SNV  1056 1047  9 112     0.991  0.903 0.945
#      segmentation      eval   SNV  1021 1015  6 144     0.994  0.876 0.931
#  no_contamination      eval   SNV  1115 1081 34  78     0.970  0.933 0.951
#            as_run      eval indel    52   41 11   7     0.788  0.854 0.820
#      segmentation      eval indel    50   40 10   8     0.800  0.833 0.816
#  no_contamination      eval indel    58   42 16   6     0.724  0.875 0.792
#            as_run eval_w100   SNV  1427 1402 25 264     0.982  0.842 0.907
#      segmentation eval_w100   SNV  1378 1363 15 303     0.989  0.818 0.896
#  no_contamination eval_w100   SNV  1485 1437 48 229     0.968  0.863 0.912
#            as_run eval_w100 indel    67   54 13  18     0.806  0.750 0.777
#      segmentation eval_w100 indel    64   53 11  19     0.828  0.736 0.779
#  no_contamination eval_w100 indel    77   57 20  15     0.740  0.792 0.765
#
# False positives only without the contamination table: 31 (26 SNVs); median tumour AF 0.042
#   gnomAD AF >= 1%: 2; rarer: 8; absent from gnomAD: 21
####
