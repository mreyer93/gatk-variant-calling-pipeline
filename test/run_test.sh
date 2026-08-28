#!/usr/bin/env bash
# End-to-end smoke test for the somatic pipeline on a small public dataset.
#
# Downloads ~2 MB from nf-core/test-datasets (sarek branch): a tiny human reference
# subset and a matched tumour/normal FASTQ pair, aligns them, then runs the full somatic
# workflow: read groups -> BQSR -> Mutect2 (tumour-vs-reference and tumour-vs-normal)
# -> orientation-bias and contamination filtering -> depth -> report.
#
# Annotation (CADD, Funcotator) is switched off: those databases are hundreds of
# gigabytes and are not needed to demonstrate calling. See manual/somatic.md.
#
# Usage:
#   ./test/run_test.sh              # download if needed, then run
#   ./test/run_test.sh --dry-run    # build the DAG only, run nothing
#
# Runtime is roughly two minutes on a laptop.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

DRY_RUN=""
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN="-n"

# Snakemake manages the conda envs by default. Set USE_CONDA=0 if gatk/samtools/bwa/R
# are already on PATH.
CONDA_FLAG="--use-conda"
[[ "${USE_CONDA:-1}" == "0" ]] && CONDA_FLAG=""

D=test/data
REF="$D/reference"; FQ="$D/fastq"; BAM="$D/bam"
RAW="https://raw.githubusercontent.com/nf-core/test-datasets/sarek"
mkdir -p "$REF" "$FQ" "$BAM"

fetch() { [[ -s "$2" ]] && return 0; echo "  $(basename "$2")"; curl -sfL -o "$2" "$1" || { echo "FAILED: $1" >&2; exit 1; }; }

echo "==> Reference"
fetch "$RAW/reference/human_g1k_v37_decoy.small.fasta" "$REF/human_g1k_v37_decoy.small.fasta"
fetch "$RAW/reference/human_g1k_v37_decoy.small.fasta.fai" "$REF/human_g1k_v37_decoy.small.fasta.fai"
fetch "$RAW/reference/dbsnp_138.b37.small.vcf.gz" "$REF/dbsnp_138.b37.small.vcf.gz"
fetch "$RAW/reference/dbsnp_138.b37.small.vcf.gz.tbi" "$REF/dbsnp_138.b37.small.vcf.gz.tbi"
fetch "$RAW/reference/gnomAD.r2.1.1.GRCh37.small.PASS.AC.AF.only.vcf.gz" "$REF/gnomad.small.vcf.gz"
fetch "$RAW/reference/gnomAD.r2.1.1.GRCh37.small.PASS.AC.AF.only.vcf.gz.tbi" "$REF/gnomad.small.vcf.gz.tbi"

# The pipeline normalises BAM contigs to chr1, chr2, ... so the reference and the
# resource VCFs have to use the same naming.
if [[ ! -s "$REF/genome.chr.fasta.fai" ]]; then
    echo "==> Renaming reference contigs to chr-prefixed"
    awk '/^>/{sub(/^>/,">chr"); print; next} {print}' \
        "$REF/human_g1k_v37_decoy.small.fasta" > "$REF/genome.chr.fasta"
    samtools faidx "$REF/genome.chr.fasta"
    rm -f "$REF/genome.chr.dict"
    gatk CreateSequenceDictionary -R "$REF/genome.chr.fasta" -O "$REF/genome.chr.dict" > /dev/null
    cut -f1 "$REF/human_g1k_v37_decoy.small.fasta.fai" | awk '{print $1"\tchr"$1}' > "$REF/rename_chrs.txt"
    bcftools annotate --rename-chrs "$REF/rename_chrs.txt" "$REF/dbsnp_138.b37.small.vcf.gz" \
        -Oz -o "$REF/dbsnp_138.chr.vcf.gz" && tabix -f -p vcf "$REF/dbsnp_138.chr.vcf.gz"
    bcftools annotate --rename-chrs "$REF/rename_chrs.txt" "$REF/gnomad.small.vcf.gz" \
        -Oz -o "$REF/gnomAD.chr.vcf.gz" && tabix -f -p vcf "$REF/gnomAD.chr.vcf.gz"
fi

echo "==> FASTQ (matched tumour/normal, multiple lanes)"
python3 - "$FQ" <<'PY'
import json, os, sys, urllib.request
fq = sys.argv[1]
base = "https://raw.githubusercontent.com/nf-core/test-datasets/sarek/testdata/tiny"
api = "https://api.github.com/repos/nf-core/test-datasets/contents/testdata/tiny/{}?ref=sarek"
for d in ("normal", "tumor"):
    for f in json.load(urllib.request.urlopen(api.format(d))):
        n = f["name"]
        if not n.endswith(".fastq.gz"):
            continue
        dest = os.path.join(fq, n)
        if os.path.exists(dest) and os.path.getsize(dest):
            continue
        print("  " + n)
        urllib.request.urlretrieve(f"{base}/{d}/{n}", dest)
PY

if [[ ! -s "$BAM/TUMOUR.bam" || ! -s "$BAM/NORMAL.bam" ]]; then
    echo "==> Aligning (bwa mem), merging lanes per sample"
    [[ -s "$REF/genome.chr.fasta.bwt" ]] || bwa index "$REF/genome.chr.fasta" 2>/dev/null
    for s in n t; do
        name=$([ "$s" = n ] && echo NORMAL || echo TUMOUR)
        parts=()
        for r1 in "$FQ"/tiny_${s}_L*_R1_xxx.fastq.gz; do
            [[ -e "$r1" ]] || continue
            r2="${r1/_R1_/_R2_}"; lane=$(basename "$r1" | sed -E 's/.*_(L[0-9]+)_.*/\1/')
            bwa mem -t "${CORES:-4}" \
                -R "@RG\tID:${lane}\tSM:${name}\tLB:lib1\tPL:ILLUMINA\tPU:${lane}" \
                "$REF/genome.chr.fasta" "$r1" "$r2" 2>/dev/null \
                | samtools sort -o "$BAM/${name}_${lane}.bam" -
            parts+=("$BAM/${name}_${lane}.bam")
        done
        samtools merge -f -o "$BAM/${name}.bam" "${parts[@]}"
        samtools index "$BAM/${name}.bam"
        rm -f "${parts[@]}"
        echo "  ${name}: $(samtools view -c "$BAM/${name}.bam") reads"
    done
fi

echo "==> Writing sample sheet and config"
printf 'sample\tbamfile\tpatient\ttimepoint\ttumor_normal\tlibrary\tplatform\tplatform_unit\n' > "$D/bam_metadata.tsv"
printf 'TUMOUR\t%s/TUMOUR.bam\tPATIENT_A\tDx\tT\tlib1\tILLUMINA\tunit1\n' "$BAM" >> "$D/bam_metadata.tsv"
printf 'NORMAL\t%s/NORMAL.bam\tPATIENT_A\tDx\tN\tlib1\tILLUMINA\tunit1\n' "$BAM" >> "$D/bam_metadata.tsv"

awk 'BEGIN{OFS="\t"} $1=="chr1"{print $1,0,$2,"+",$1"_region"}' "$REF/genome.chr.fasta.fai" > "$REF/targets.bed"

cat > "$D/config_test.yaml" <<EOF
project_name: "Somatic calling demo (nf-core/sarek tiny tumour-normal pair)"
output_directory: "$D/results"
bam_metadata: "$D/bam_metadata.tsv"
annotation_only: False
skip_bqsr: False
skip_annotation: True
call_tumor_normal: True
filter_min_af: 0.02
filter_min_dp: 5
filter_GATK: ['PASS', 'normal_artifact']
filter_min_cadd: 20
target_bed_w100: "$REF/targets.bed"
genome_version: "GRCh37"
REF_FILE: "$REF/genome.chr.fasta"
funcotator_data_path: "not_used"
dbsnp_file: "$REF/dbsnp_138.chr.vcf.gz"
dbsnp_common_file: "$REF/dbsnp_138.chr.vcf.gz"
gnomad_file: "$REF/gnomAD.chr.vcf.gz"
make_report: True
report_pdf: False
EOF

echo "==> Running pipeline"
snakemake -s call_bam_GATK/call_bam_GATK.snakefile \
    --configfile "$D/config_test.yaml" --cores "${CORES:-4}" $CONDA_FLAG $DRY_RUN

if [[ -z "$DRY_RUN" ]]; then
    echo; echo "==> Checking expected outputs"
    fail=0
    for f in "$D/results/02_variants_reference/01_gatk_variant_calling/filtered/TUMOUR_filtered.vcf" \
             "$D/results/03_variants_TvN/01_gatk_variant_calling/filtered/PATIENT_A_Dx_filtered.vcf" \
             "$D/results/01_prepare_bam/depth_smoothed_join.tsv" \
             "$D/results/09_report/somatic_report.html"; do
        if [[ -s "$f" ]]; then echo "  OK   $f"; else echo "  MISS $f"; fail=1; fi
    done
    echo
    if [[ $fail -eq 0 ]]; then
        t=$(awk '!/^#/ && $7=="PASS"' "$D/results/02_variants_reference/01_gatk_variant_calling/filtered/TUMOUR_filtered.vcf" | wc -l | tr -d ' ')
        p=$(awk '!/^#/ && $7=="PASS"' "$D/results/03_variants_TvN/01_gatk_variant_calling/filtered/PATIENT_A_Dx_filtered.vcf" | wc -l | tr -d ' ')
        echo "Smoke test PASSED - $t PASS calls tumour-vs-reference, $p after subtracting the matched normal."
    else
        echo "Smoke test FAILED - see $D/results/logs/" >&2; exit 1
    fi
fi
