gatk FilterMutectCalls \
    -R /path/to/data/references/human_g1k_v37_fixedChr.fasta  \
    -V output/03_variants_TvN/001-103_unfiltered.vcf \
    -stats output/03_variants_TvN/001-103_unfiltered.vcf.stats \
    -O output/03_variants_TvN/001-103_filtered.vcf

gatk Funcotator \
    --variant output/03_variants_TvN/001-103_filtered.vcf \
    --reference /path/to/data/references/human_g1k_v37_fixedChr.fasta \
    --ref-version hg19 \
    --data-sources-path /path/to/data/references/funcotator_dataSources.v1.7.20200521s/ \
    --output output/funcotator_TvN/001-103_filtered.vcf \
    --output-file-format VCF

for i in *; do echo $i; j=$(echo $i | sed "s/DS-[0-9][0-9][0-9][0-9][0-9][0-9]_//"); mv $i $j; done


cat GRCh37_af-only-gnomad.raw.sites.vcf | \
    sed -e 's/^1/chr1/' | sed -e 's/^2/chr2/' | \
    sed -e 's/^3/chr3/' | sed -e 's/^4/chr4/' | \
    sed -e 's/^5/chr5/' | sed -e 's/^6/chr6/' | \
    sed -e 's/^7/chr7/' | sed -e 's/^8/chr8/' | \
    sed -e 's/^9/chr9/' | sed -e 's/^10/chr10/' | \
    sed -e 's/^11/chr11/' | sed -e 's/^12/chr12/' | \
    sed -e 's/^13/chr13/' | sed -e 's/^14/chr14/' | \
    sed -e 's/^15/chr15/' | sed -e 's/^16/chr16/' | \
    sed -e 's/^17/chr17/' | sed -e 's/^18/chr18/' | \
    sed -e 's/^19/chr19/' | sed -e 's/^20/chr20/' | \
    sed -e 's/^21/chr21/' | sed -e 's/^22/chr22/' | \
    sed -e 's/^X/chrX/' | sed -e 's/^Y/chrY/' | \
    sed -e 's/^MT/chrM/' 

003-101_Dx_BM
003-101_Dx_PB
003-101_Dx_BU

gatk Mutect2 \
-R /path/to/data/references/human_g1k_v37_fixedChr.fasta \
-I output/01_prepare_bam/recalibrate/003-101_Dx_BM.fixchr.bam \
-I output/01_prepare_bam/recalibrate/003-101_Dx_PB.fixchr.bam \
-I output/01_prepare_bam/recalibrate/003-101_Dx_BU.fixchr.bam \
-normal 003-101_Dx_BU \
-O output/03_variants_TvN/unfiltered/003-101_Dx_unfiltered.vcf \
-L /path/to/data/references/AML_261_IDT/targets_w100_hg19.st.bed \
--native-pair-hmm-threads 4 \
--f1r2-tar-gz output/03_variants_TvN/unfiltered/003-101_Dx_f1r2.tar.gz
    

BAM=output/01_prepare_bam/recalibrate/008-103_Dx_BU.fixchr.bam
ls *fixchr.bam | xargs -P 16 -I {} sh -c "echo {}; samtools view -H {} | sed "s/DS-[0-9][0-9][0-9][0-9][0-9][0-9]_//" | samtools reheader - {} > ../newhead/{} && mv ../newhead/{} {}"

# remove DS from the f1r2.tar.gz files
samps=( 003-101_Dx 007-103_Dx 007-104_Dx 007-104_ITPD84 007-104_ITPEOW 008-102_Dx 008-102_ITPD84 008-104_Dx 008-104_EOT 009-101_Dx 010-101_Dx 010-103_Dx 011-101_Dx 012-102_Dx 014-101_Dx 020-102_Dx 060-105_Dx 060-107_Dx )
for s in "${samps[@]}"; do
    echo $s
    tar xvfz "$s"_f1r2.tar.gz
    for i in DS*"$s"*; do echo $i; sed -i "s/DS-[0-9][0-9][0-9][0-9][0-9][0-9]_//" $i ; j=$(echo $i | sed "s/DS-[0-9][0-9][0-9][0-9][0-9][0-9]_//"); mv $i $j; done
    tar -czvf "$s"_f1r2.tar.gz "$s"*.alt* "$s"*.ref*
    rm "$s"*.alt* "$s"*.ref*
done


for i in *.vcf; do echo $i; sed -i "s/DS-[0-9][0-9][0-9][0-9][0-9][0-9]_//" $i; done
for i in *.vcf; do echo $i;  gatk IndexFeatureFile -I $i;done

DS-358116_008-103_Dx_BS

        gatk Mutect2             -R /path/to/data/references/human_g1k_v37_fixedChr.fasta             -I output/01_prepare_bam/recalibrate/011-108_Dx_BS.fixchr.bam -I output/01_prepare_bam/recalibrate/011-108_ITPD84_PG.fixchr.bam             -normal 011-108_Dx_BS             -O output/03_variants_TvN/unfiltered/011-108_ITPD84_unfiltered.vcf             -L /path/to/data/references/AML_261_IDT/targets_w100_hg19.st.bed             --native-pair-hmm-threads 4             --f1r2-tar-gz output/03_variants_TvN/unfiltered/011-108_ITPD84_f1r2.tar.gz
   

# sort by last column
 f=007-101_Dx_filtered_annotated.vcf; awk '{print $NF,$0}' "$f" | sort -nr | cut -f2- -d' ' | catcol | less -S


 # redo bam headers 
 A USER ERROR has occurred: Bad input: Sample 008-103_Dx_BS is not in BAM header: [008-103_ATP3D168_PG, 008-103_Dx_BU, DS-358116_008-103_Dx_BS]                                                                                                                                           
