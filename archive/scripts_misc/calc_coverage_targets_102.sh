# basically the same as the processing snakefile but done as a loop
# over many files. Here for the ega_data
samp=$1
echo "STARTING $samp"
ref=~/data/references/human_g1k_v37_fixedChr.fasta.fai
IDT_target=~/data/references/AML_261_IDT/targets_w100_hg19.st.bed
rl=75
in_bam=/path/to/data/ctp_102_calling/gatk_calling_bam/output/01_prepare_bam/recalibrate/"$samp".fixchr.bam
out_cov=/path/to/data/ctp_102_calling/gatk_calling_bam/output/01_prepare_bam/cov/"$samp".w100.cov
out_cov_stats=/path/to/data/ctp_102_calling/gatk_calling_bam/output/01_prepare_bam/cov/"$samp".w100.cov.stats

cut -f 1,2 "$ref" > "$out_cov".genome
coverageBed -sorted -g "$out_cov".genome -a "$IDT_target" -b "$in_bam" > "$out_cov".tmp
awk -v r="$rl"  '{ print $6*r/$8; }' "$out_cov".tmp | paste "$out_cov".tmp  - > "$out_cov"
echo "$out_cov"  > "$out_cov_stats"
stats -c 8 "$out_cov"  >> "$out_cov_stats"
rm "$out_cov".tmp "$out_cov".genome
