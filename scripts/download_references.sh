#!/usr/bin/env bash
# Downloads the GRCh38 reference resources used by call_bam_GATK/ and process_reads/
# from Broad Institute's public GCS buckets (no auth needed for these specific buckets).
#
# This does NOT cover:
#   - Funcotator data sources: run `gatk FuncotatorDataSourceDownloader` (see bottom of
#     this script's output) - the hosting location changed in GATK 4.6.2, so a static
#     URL here would go stale.
#   - CADD's scored-variant databases (400G+): run `cadd-install.sh` from the
#     `gatk-pipeline` conda env (see manual/requirements.md) - it's an interactive
#     installer, not a simple file fetch.
#   - Your own capture-panel target BED file - there is no generic default.
#
# Usage:
#   ./scripts/download_references.sh [output_dir]      # download everything (~17G)
#   ./scripts/download_references.sh --dry-run [output_dir]   # show what would happen, no download
#
# Safe to re-run: skips files that already exist and look complete (checked by size).
# Downloads are resumable (curl -C -) since some of these files are multi-gigabyte.

set -euo pipefail

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    shift
fi
OUTDIR="${1:-references}"

BASE_REF="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
BASE_SOMATIC="https://storage.googleapis.com/gatk-best-practices/somatic-hg38"

# name, url, approx size (for the dry-run summary only - not enforced)
FILES=(
    "Homo_sapiens_assembly38.fasta|$BASE_REF/Homo_sapiens_assembly38.fasta|~3.1G"
    "Homo_sapiens_assembly38.fasta.fai|$BASE_REF/Homo_sapiens_assembly38.fasta.fai|~20K"
    "Homo_sapiens_assembly38.dict|$BASE_REF/Homo_sapiens_assembly38.dict|~10K"
    "dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf|$BASE_REF/Homo_sapiens_assembly38.dbsnp138.vcf|~10.4G"
    "dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf.idx|$BASE_REF/Homo_sapiens_assembly38.dbsnp138.vcf.idx|~2.5M"
    "af-only-gnomad.hg38.vcf.gz|$BASE_SOMATIC/af-only-gnomad.hg38.vcf.gz|~3.0G"
    "af-only-gnomad.hg38.vcf.gz.tbi|$BASE_SOMATIC/af-only-gnomad.hg38.vcf.gz.tbi|~2M"
    "dbsnp/small_exac_common_3.hg38.vcf.gz|$BASE_SOMATIC/small_exac_common_3.hg38.vcf.gz|~1M"
    "dbsnp/small_exac_common_3.hg38.vcf.gz.tbi|$BASE_SOMATIC/small_exac_common_3.hg38.vcf.gz.tbi|~30K"
)

echo "Reference download plan (output dir: $OUTDIR)"
echo "-----------------------------------------------"
for entry in "${FILES[@]}"; do
    IFS='|' read -r name url size <<< "$entry"
    printf '  %-55s %s\n' "$name" "$size"
done
echo "-----------------------------------------------"
echo "Total: ~17G. This maps onto config_call_bam_GATK.yaml / config_processing.yaml as:"
echo "  REF_FILE / reference_file  -> $OUTDIR/Homo_sapiens_assembly38.fasta"
echo "  dbsnp_file                 -> $OUTDIR/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf"
echo "  gnomad_file                -> $OUTDIR/af-only-gnomad.hg38.vcf.gz"
echo "  dbsnp_common_file          -> $OUTDIR/dbsnp/small_exac_common_3.hg38.vcf.gz"
echo ""

if $DRY_RUN; then
    echo "(dry run - nothing downloaded)"
    exit 0
fi

mkdir -p "$OUTDIR/dbsnp"

for entry in "${FILES[@]}"; do
    IFS='|' read -r name url size <<< "$entry"
    dest="$OUTDIR/$name"
    if [[ -s "$dest" ]]; then
        echo "Already present, skipping: $dest"
        continue
    fi
    echo "Downloading $name ($size)..."
    curl -L -C - --fail -o "$dest" "$url"
done

echo ""
echo "Done. Remaining manual steps:"
echo "  1. Funcotator data sources (GRCh38) - run inside the gatk-pipeline conda env:"
echo "       gatk FuncotatorDataSourceDownloader --germline --hg38 --validate-integrity --extract-after-download"
echo "       gatk FuncotatorDataSourceDownloader --somatic  --hg38 --validate-integrity --extract-after-download"
echo "  2. CADD databases (400G+) - run inside the gatk-pipeline conda env: cadd-install.sh"
echo "  3. Your capture-panel target BED file(s) - no generic default, bring your own per project."
