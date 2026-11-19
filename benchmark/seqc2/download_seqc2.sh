#!/usr/bin/env bash
# Downloads the SEQC2 HCC1395 benchmark data used by manual/benchmark.md, about 27 GB:
#   - the Illumina-centre whole-exome BAMs (tumour WES_IL_T_1, 12.5 GB; normal WES_IL_N_1,
#     14.1 GB), BWA-MEM alignments to GRCh38 with duplicates marked, checked against the
#     published md5 sums
#   - the Agilent SureSelect V6 + UTR capture targets lifted to hg38, and the padded calling
#     region the benchmark run uses (targets +/- 100 bp, primary chromosomes)
#   - the SEQC2 somatic reference call set v1.2.1, high-confidence regions and the all-caller
#     supersets (truth/)
# Source: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/
# (Fang et al. 2021, Nat Biotechnol 39:1151; Xiao et al. 2021, Nat Biotechnol 39:1141).
# Re-running skips files already downloaded.
#
# Usage: benchmark/seqc2/download_seqc2.sh <destination dir>

set -euo pipefail
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <destination dir>" >&2
    exit 1
fi
dest=$1
base=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG
mkdir -p "$dest/truth"
cd "$dest"

get() {  # url, output file
    if [[ -s "$2" ]]; then echo "    have $2"; return; fi
    echo "    fetching $2"
    curl -fL --retry 5 --retry-delay 10 -C - -o "$2.part" "$1" && mv "$2.part" "$2"
}

echo "==> Exome BAMs"
for s in WES_IL_T_1 WES_IL_N_1; do
    get "$base/data/WES/$s.bwa.dedup.bam" "$s.bwa.dedup.bam"
    get "$base/data/WES/$s.bwa.dedup.bai" "$s.bwa.dedup.bai"
done
get "$base/data/WES/WES.md5" WES.md5
echo "==> Checking md5 sums"
for s in WES_IL_T_1 WES_IL_N_1; do
    want=$(awk -v f="WES/alignment/$s.bwa.dedup.bam" '$2 == f {print $1}' WES.md5)
    have=$( (command -v md5sum > /dev/null && md5sum "$s.bwa.dedup.bam" || md5 -r "$s.bwa.dedup.bam") | cut -d' ' -f1)
    if [[ "$want" != "$have" ]]; then
        echo "md5 mismatch for $s.bwa.dedup.bam (expected $want, got $have)" >&2
        exit 1
    fi
    echo "    $s ok"
done

echo "==> Capture targets"
bed=S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed6.gz
get "$base/technical/reference_genome/Exome_Target_bed/$bed" "$bed"
# Padded calling region for the pipeline's target_bed_w100: primary chromosomes, 100 bp
# either side, gene name in column 5, sorted by chromosome and position
gzip -dc "$bed" | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ {s=$2-100; if (s<0) s=0; split($4, n, ","); print $1, s, $3+100, $6, n[1]}' \
    | sort -s -k1,1V -k2,2n -k3,3n > V6_UTR_w100.bed

echo "==> Truth set"
for f in README.md High-Confidence_Regions_v1.2.bed \
         high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz \
         sSNV.MSDUKT.superSet.v1.2.vcf.gz sINDEL.MDKT.superSet.v1.2.vcf.gz; do
    get "$base/release/latest/$f" "truth/$f"
done
for f in truth/*.vcf.gz; do
    [[ -s "$f.tbi" ]] || tabix -p vcf "$f"
done
echo "Done: $dest"
