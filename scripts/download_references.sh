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
#   ./scripts/download_references.sh [output_dir]                        # core refs (~17G)
#   ./scripts/download_references.sh --with-alphamissense [output_dir]   # core refs + AlphaMissense (~17.6G)
#   ./scripts/download_references.sh --dry-run [output_dir]              # show what would happen, no download
#
# AlphaMissense (germline missense deleteriousness, ~640M) is opt-in since it's only
# useful for the germline pipeline and isn't needed for a CADD-only setup. Requires
# `tabix` on PATH (from the gatk-pipeline conda env's htslib package) to index it.
#
# Safe to re-run: skips files that already exist and look complete (checked by size).
# Downloads are resumable (curl -C -) since some of these files are multi-gigabyte.

set -euo pipefail

DRY_RUN=false
WITH_AM=false
ARGS=()
for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=true ;;
        --with-alphamissense) WITH_AM=true ;;
        *) ARGS+=("$arg") ;;
    esac
done
OUTDIR="${ARGS[0]:-references}"

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

if $WITH_AM; then
    FILES+=("AlphaMissense_hg38.tsv.gz|https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz?download=1|~640M")
fi

echo "Reference download plan (output dir: $OUTDIR)"
echo "-----------------------------------------------"
for entry in "${FILES[@]}"; do
    IFS='|' read -r name url size <<< "$entry"
    printf '  %-55s %s\n' "$name" "$size"
done
echo "-----------------------------------------------"
TOTAL_SIZE="~17G"
if $WITH_AM; then TOTAL_SIZE="~17.6G"; fi
echo "Total: $TOTAL_SIZE. This maps onto config_call_bam_GATK.yaml / config_processing.yaml as:"
echo "  REF_FILE / reference_file  -> $OUTDIR/Homo_sapiens_assembly38.fasta"
echo "  dbsnp_file                 -> $OUTDIR/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf"
echo "  gnomad_file                -> $OUTDIR/af-only-gnomad.hg38.vcf.gz"
echo "  dbsnp_common_file          -> $OUTDIR/dbsnp/small_exac_common_3.hg38.vcf.gz"
if $WITH_AM; then
    echo "  alphamissense_file (germline config) -> $OUTDIR/AlphaMissense_hg38.tsv.gz"
fi
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

if $WITH_AM; then
    am_file="$OUTDIR/AlphaMissense_hg38.tsv.gz"
    if [[ -s "${am_file}.tbi" ]]; then
        echo "Already indexed, skipping: ${am_file}.tbi"
    else
        echo "Indexing AlphaMissense scores (tabix)..."
        if ! command -v tabix &> /dev/null; then
            echo "ERROR: tabix not found - activate the gatk-pipeline conda env first (it ships with htslib)." >&2
            exit 1
        fi
        # already BGZF-compressed (confirmed via its gzip header), so tabix can index it directly
        tabix -s 1 -b 2 -e 2 -c '#' "$am_file"
    fi
fi

echo ""
echo "Done. Remaining manual steps:"
echo "  1. Funcotator data sources (GRCh38) - run inside the gatk-pipeline conda env:"
echo "       gatk FuncotatorDataSourceDownloader --germline --hg38 --validate-integrity --extract-after-download"
echo "       gatk FuncotatorDataSourceDownloader --somatic  --hg38 --validate-integrity --extract-after-download"
echo "  2. CADD databases (400G+) - run inside the gatk-pipeline conda env: cadd-install.sh"
echo "  3. Your capture-panel target BED file(s) - no generic default, bring your own per project."
