s=PDNEW36a
ref=/path/to/data/references/human_g1k_v37_fixedChr.fasta
intervals=/path/to/data/references/sure_select_V5/S04380110_Padded.bed
inbam=/path/to/data/EGA_data/01_prepare_bam/recalibrate/"$s".fixchr.bam
outgvcf=/path/to/data/EGA_data/05_haplotypecaller/"$s".g.vcf.gz
gatk HaplotypeCaller --java-options "-Xmx7g" \
    -R "$ref" \
    -I "$inbam" \
    -O "$outgvcf" \
    -ERC GVCF \
    -L "$intervals" \
    --native-pair-hmm-threads 2 

# NEED TO RENAME THE SAMPLE NAMES IN THE VCFS WITH TAHT SCRIPT AGAIN

ref=/path/to/data/references/human_g1k_v37_fixedChr.fasta
intervals=/path/to/data/references/sure_select_V5/S04380110_Padded.bed


# what if I try iterative import into the one that exists already? 
gatk GenomicsDBImport \
    --reference "$ref" \
    -V 1711STDY5270146.g.vcf.gz \
    --genomicsdb-update-workspace-path PD477DB \
    --reader-threads 4 \
    --max-num-intervals-to-import-in-parallel 4 \
    -L chr1

gatk --java-options "-Xmx200g" GenomicsDBImport \
    --reference "$ref" \
    --sample-name-map sample_list.txt \
    --genomicsdb-workspace-path ../06_GDB/chr1 \
    --reader-threads 32 \
    -L chr1

gatk --java-options "-Xmx200g" GenomicsDBImport \
    --reference "$ref" \
    --sample-name-map sample_list.txt \
    --genomicsdb-workspace-path ../06_GDB/chr2 \
    --reader-threads 16 \
    -L chr2

gatk --java-options "-Xmx20g" GenomicsDBImport \
    --reference "$ref" \
    --sample-name-map sample_list.txt \
    --genomicsdb-workspace-path ../06_GDB/chrX \
    --reader-threads 16 \
    -L chrX

gatk --java-options "-Xmx10g" GenomicsDBImport \
    --reference /path/to/data/references/human_g1k_v37_fixedChr.fasta \
    --sample-name-map sample_list.txt \
    --genomicsdb-workspace-path ../06_GDB/chr{} \
    --reader-threads 2 \
    --batch-size 32 \
    -L chr{}

gatk GenotypeGVCFs \
    --reference /path/to/data/references/human_g1k_v37_fixedChr.fasta \
    --variant gendb://../06_GDB/chr1 \
    --output ../06_GDB/chr1.vcf 


gatk --java-options "-Xmx10g" GenotypeGVCFs \
    --reference /path/to/data/references/human_g1k_v37_fixedChr.fasta \
    --variant gendb://../06_GDB/chr{} \
    --output ../06_GDB/chr{}.vcf 

gatk SelectVariants \
    --reference /path/to/data/references/human_g1k_v37_fixedChr.fasta \
    --variant chr{}.vcf \
    --select-type-to-include SNP \
    --output chr{}_selected.vcf \
    -L ~/data/references/sure_select_V5/S04380110_Padded.bed

# variant recalibrator to look for high quiality variants
echo "$(seq 1 22 ) X" | tr " " "\n" | xargs -I {} -P 32 sh -c "gatk VariantRecalibrator \
   -R /path/to/data/references/human_g1k_v37_fixedChr.fasta \
   -V chr{}_selected.vcf \
   -O recal/chr{}.recal \
   --tranches-file tranches/chr{}.tranches \
   --rscript-file outR/chr{}.R \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /path/to/data/references/vqsr_datasets/hapmap_3.3.b37.chr.vcf.gz \
   --resource:omni,known=false,training=true,truth=false,prior=12.0 /path/to/data/references/vqsr_datasets/1000G_omni2.5.b37.chr.vcf.gz \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 /path/to/data/references/vqsr_datasets/1000G_phase1.snps.high_confidence.b37.chr.vcf.gz \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /path/to/data/references/vqsr_datasets/dbsnp_138.b37.chr.vcf.gz \
   --ignore-all-filters true
"

# apply the quality score recalibration
echo "$(seq 1 22 ) X" | tr " " "\n" | xargs -I {} -P 32 sh -c "gatk ApplyVQSR \
   -R /path/to/data/references/human_g1k_v37_fixedChr.fasta \
   -V chr{}_selected.vcf \
   -O chr{}_selected_vqsr.vcf \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file tranches/chr{}.tranches \
   --recal-file recal/chr{}.recal \
   -mode SNP \
   --ignore-all-filters true
 "

# Funcotation and cadd score
echo "$(seq 1 22 ) X" | tr " " "\n" | xargs -I {} -P 32 sh -c "gatk Funcotator \
    --variant chr{}_selected_vqsr.vcf \
    --reference /path/to/data/references/human_g1k_v37_fixedChr.fasta \
    --ref-version hg19 \
    --data-sources-path /path/to/data/references/funcotator_dataSources.v1.7.20200521s/ \
    --output chr{}_selected_vqsr_func.vcf \
    --output-file-format VCF
    "


# CADD scoring 
sed "s/^chr//g" chr{}_selected_filtered.vcf > /tmp/chr{}.vcf && \
CADD.sh -a -g GRCh37 -o chr{}_CADD_full.tsv.gz /tmp/chr{}.vcf && \
zcat -f chr{}_CADD_full.tsv | cut -f 1,2,3,4,5,6,10,21,116 > chr{}_CADD_simple.tsv && \
rm /tmp/chr{}.vcf 


# combine annotations with scripts/combine_annotation.R ??? 
mkdir decomposed
echo "$(seq 1 22 ) X" | tr " " "\n" | xargs -I {} -P 32 sh -c "vt decompose chr{}_selected_vqsr_func.vcf > decomposed/chr{}_decomposed.vcf && \
  vcftools --vcf decomposed/chr{}_decomposed.vcf --extract-FORMAT-info GT --out decomposed/chr{}_decomposed && \
  vcftools --vcf decomposed/chr{}_decomposed.vcf --extract-FORMAT-info DP --out decomposed/chr{}_decomposed
  "

mkdir annotation_combined
echo "$(seq 1 22 ) X" | tr " " "\n" | xargs -I {} -P 32 sh -c "Rscript ~/projects/project/EGA_variant_processing/germline/combine_annotation_germline.R \
  decomposed/chr{}_decomposed.vcf \
  decomposed/chr{}_decomposed.GT.FORMAT \
  chr{}_CADD_simple.tsv \
  ~/data/EGA_data/bam_metadata_uniq_samp.tsv \
  annotation_combined/chr{}_combined_unfiltered.vcf \
  annotation_combined/chr{}_combined_filtered.vcf
"




###################################################################################################################
# # # OLD CODE ####################################################################################################
###################################################################################################################
tail -n +2 sample_list.txt | head -n -1 | xargs -I {} -P 20 sh -c "./xargs_haplotype_script.sh {}"
cat ../../sample_list.txt | xargs -I {} -P 16 sh -c "echo {}; gatk SelectVariants -V {}_filtered.vcf.gz -L ~/data/references/sure_select_V5/S04380110_Padded.bed -O ../filtered_selected/{}_filtered.vcf.gz"
cat ../../sample_list.txt | xargs -I {} -P 16 sh -c "echo {}; gatk SelectVariants -V {}_filtered_fucnc.vcf.gz -L ~/data/references/sure_select_V5/S04380110_Padded.bed -O ../annotate_Funcotator_selected/{}_filtered_func.vcf.gz"

# depth of 5 in at least one normal sample
# remove key column 
# rounding on score new
BOX passwords
vNXQj53TiDFgaQU

# HARDFILT manual filtration that I no longer do
# gatk VariantFiltration  \
#     --reference /path/to/data/references/human_g1k_v37_fixedChr.fasta \
#     --variant chr{}_selected.vcf \
#     --output chr{}_selected_filtered.vcf \
#     --filter-expression "QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
#     --filter-name "gatk_hardFilt"
