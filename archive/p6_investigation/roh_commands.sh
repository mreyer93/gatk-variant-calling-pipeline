###################################################################################################################################
############################# 2022-07-10 ##########################################################################################
###################################################################################################################################
# change AS_FilterStatus,Number=A to AS_FilterStatus,Number=. 
sed "s/AS_FilterStatus,Number=A/AS_FilterStatus,Number=./g"

# Do for all vcf files
ls *.vcf.gz | xargs -I {} -P 8 sh -c "echo {}; zcat {} | sed 's/AS_FilterStatus,Number=A/AS_FilterStatus,Number=./g' | bgzip > {}.new; mv {}.new {}; tabix {}"
bcftools merge *.vcf.gz -Oz --threads 8 > merged.vcf.gz

# roh test on all? 
ls *.vcf.gz |  xargs -I  {} -P 8 sh -c "echo {}; bcftools roh -G30 --AF-dflt 0.4 {} > {}.roh"
ls *.roh |  xargs -I  {} -P 1 sh -c "echo {}; grep -v '#' {}|  cut -f 5| uniq "

# remove chr from all
# nochr.awk is 
# {gsub(/^chr/,""); print}
ls *.vcf.gz | xargs  -I {} -P 8 sh -c "echo {}; zcat {} | awk -f nochr.awk | bgzip >  fixchr/{}; tabix fixchr/{}"

# prepare working version of 1kg AF list
zcat AFs.tab.gz | cut -f 1,2 > col12.txt
zcat AFs.tab.gz | cut -f 3 | cut -f 1 -d "," > col3.txt
zcat AFs.tab.gz | cut -f 3 | cut -f 2-4 -d "," > col4.txt
zcat AFs.tab.gz | cut -f 4 > col5.txt
paste col12.txt col3.txt col4.txt col5.txt | bgzip > AFs_MOD.tab.gz
tabix -s1 -b2 -e2 AFs_MOD.tab.gz

bcftools annotate -c CHROM,POS,REF,ALT,AF1KG -h ../1000GP-AFs/AFs.tab.gz.hdr -a ../1000GP-AFs/AFs_MOD.tab.gz PD4060a_filtered.vcf.gz | \
   bcftools roh --AF-tag AF1KG -M 100 -m ../genetic-map/genetic_map_chr{CHROM}_combined_b37.txt -o PD4060a_filtered_roh.txt -G30

ls *.vcf.gz | xargs  -I {} -P 8 sh -c "echo {}; bcftools annotate -c CHROM,POS,REF,ALT,AF1KG -h ../1000GP-AFs/AFs.tab.gz.hdr -a ../1000GP-AFs/AFs_MOD.tab.gz {} | bcftools roh --AF-tag AF1KG -M 100 -m ../genetic-map/genetic_map_chr{CHROM}_combined_b37.txt -o {}.roh -G30"

bcftools annotate --rename-chrs test/chr-names.txt './PD7281a_filtered.vcf.gz' -Ou | \
  bcftools annotate -c CHROM,POS,REF,ALT,AF1KG -h 1000GP-AFs/AFs.tab.gz.hdr -a 1000GP-AFs/AFs.tab.gz -Ob -o test/PD7281a_filtered.bcf.part && mv test/PD7281a_filtered.bcf.part test/PD7281a_filtered.bcf

bcftools roh --AF-tag AF1KG -M 100 -m genetic-map/genetic_map_chr{CHROM}_combined_b37.txt -o roh.txt -G30 test/PD6548b_filtered.bcf

###################################################################################################################################
############################# 2022-07-11 ##########################################################################################
###################################################################################################################################
# Trying with some 1KG data
# look to see if these features are visible in the samples with bcftools roh

bcftools query -l ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz > sample_names.txt
cat sample_names.txt | xargs -I {} -P 8 sh -c "echo  {} ; bcftools roh --samples {} --threads 1 ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -G30 > roh/{}_roh.txt"
bcftools roh --samples HG00106,HG00107 --threads 8 ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -G30 > test_roh.txt

###################################################################################################################################
############################# 2022-07-12 ##########################################################################################
###################################################################################################################################
# tyring to find the "p6" pattern in the EGA data 
bcftools view --regions-file p6_pattern.bed merged.vcf.gz   > merged_p6_pattern.vcf
plink2 --vcf merged_p6_pattern.vcf --make-bfile --out merged_p6_pattern
cat merged_p6_pattern.vcf | awk '{gsub(/^chr/,""); print}'  > merged_p6_pattern_fixchr.vcf
bgizip merged_p6_pattern_fixchr.vcf
tabix merged_p6_pattern_fixchr.vcf.gz
bcftools annotate -c CHROM,POS,REF,ALT,AF1KG -h 1000GP-AFs/AFs.tab.gz.hdr -a 1000GP-AFs/AFs_MOD.tab.gz merged_p6_pattern_fixchr.vcf > merged_p6_pattern_fixchr_annotated.vcf

# also look at GGCC pattern
mkdir EGA_haplotypecaller/p6_pattern
mkdir EGA_haplotypecaller/GGCC_pattern
bcftools view --regions-file GGCC_pattern.bed EGA_haplotypecaller/chr9/chr9.vcf.gz  > EGA_haplotypecaller/GGCC_pattern/GGCC_pattern.vcf
bcftools view --regions-file p6_pattern.bed EGA_haplotypecaller/chr9/chr9.vcf.gz  > EGA_haplotypecaller/p6_pattern/p6_pattern.vcf
plink2 --vcf EGA_haplotypecaller/GGCC_pattern/GGCC_pattern.vcf --make-bfile --out EGA_haplotypecaller/GGCC_pattern/GGCC_pattern
plink2 --vcf EGA_haplotypecaller/p6_pattern/p6_pattern.vcf --make-bfile --out EGA_haplotypecaller/p6_pattern/p6_pattern

# and for 1KG
mkdir ../1KG_data/p6_pattern
mkdir ../1KG_data/GGCC_pattern
bcftools view --regions-file GGCC_pattern_nochr.bed ../1KG_data/vcf_intersect_EGA/chr9.vcf.gz  > ../1KG_data/GGCC_pattern/GGCC_pattern.vcf
bcftools view --regions-file p6_pattern_nochr.bed ../1KG_data/vcf_intersect_EGA/chr9.vcf.gz  > ../1KG_data/p6_pattern/p6_pattern.vcf
plink2 --vcf ../1KG_data/GGCC_pattern/GGCC_pattern.vcf --make-bfile --out ../1KG_data/GGCC_pattern/GGCC_pattern
plink2 --vcf ../1KG_data/p6_pattern/p6_pattern.vcf --make-bfile --out ../1KG_data/p6_pattern/p6_pattern

# no homozygotes at these 6 positions in the EGA data???????
# maybe I should re-genotype them? 
ref=~/data/references/human_g1k_v37_fixedChr.fasta
bed=~/data/p6/p6_pattern_alts.bed
bam=~/data/EGA_data/01_prepare_bam/recalibrate/PD4060a.fixchr.bam
bcftools mpileup -d 99999 -Ov -R $bed -f $ref $bam | bcftools call -c - > test.out

# for all samples 
bam_base=~/data/EGA_data/01_prepare_bam/recalibrate
cat EGA_sample_list.txt | xargs -I {} -P 16 sh -c "echo {}; bcftools mpileup -d 99999 -Ov -R $bed -f $ref $bam_base/{}.fixchr.bam | bcftools call -c - > EGA_mpileup_out/{}.vcf; bgzip EGA_mpileup_out/{}.vcf; tabix EGA_mpileup_out/{}.vcf.gz"

# merge these all together, put into plink, load in python
cat EGA_sample_list.txt | xargs -I {} -P 16 sh -c "echo {}; bgzip EGA_mpileup_out/{}.vcf; tabix EGA_mpileup_out/{}.vcf.gz"
bcftools merge -Oz EGA_mpileup_out/*.vcf.gz > EGA_mpileup_out/merged.vcf.gz
plink2 --vcf EGA_mpileup_out/merged.vcf.gz --make-bfile --out EGA_mpileup_out/merged

# 2022-07-13 working with haplotypecaller outputs on the server
# Can we detect ROH in the haplotype called vcf? 
bcftools annotate --rename-chrs rename_chrs.txt 'chr9.vcf.gz' -Ou | \
  bcftools annotate -c CHROM,POS,REF,ALT,AF1KG -h ~/data/p6/AFs.tab.gz.hdr -a ~/data/p6/AFs_chr9.tab.gz -Oz --threads 8 -o chr9_1KG.vcf.gz 
bcftools roh --samples PD4060a --threads 8 -G30 > PD4060a_roh.txt

### annotation didn't work, use AF tag instead
bcftools query -l chr9.vcf.gz > sample_list.txt
cat sample_list.txt | xargs -I {} -P 16 sh -c "echo {}; bcftools roh --samples {} --threads 8 -G30 chr9.vcf.gz > {}_roh.txt"

bcftools roh --samples PDNEW36a --threads 8 -G30 chr9.vcf.gz > PDNEW36a_roh.txt

###################################################################################################################################
############################# 2022-07-15 ##########################################################################################
###################################################################################################################################
# Find p6 pattern in the 1KG data
# working on local
bcftools view --regions-file p6_pattern_nochr.bed ~/local_data/1KG/individual_chromosomes/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz > 1KG_p6_pattern.vcf
plink2 --vcf 1KG_p6_pattern.vcf --make-bfile --out 1KG_p6_pattern
# convert to genotype matrix with pandas_plink


# run the ROH detection for all chromosomes 
# prepare bgzipped and indexed version to use for each chr
seq 1 22 | xargs -I {} -P 8 sh -c "echo {}; mkdir chr{}; cat ~/data/EGA_data_whole_exome/06_GDB/chr{}.vcf | bgzip --threads 4 > chr{}/chr{}.vcf.gz; tabix chr{}/chr{}.vcf.gz"

for i in $(seq 1 22); do
  echo $i
  cat sample_list.txt | xargs -I {} -P 16 sh -c "echo {}; bcftools roh --samples {} --threads 2 -G30 chr$i/chr$i.vcf.gz > chr$i/{}_roh.txt"
done

# download the 1KG data for each chromosome to run the roh calculations on it. 
cd /path/to/data/1KG_data
seq 1 22 | xargs -I {} -P 8 sh -c "wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
seq 1 22 | xargs -I {} -P 8 sh -c "wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"

# run the ROH calculation for all of it 
bcftools query -l vcf/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz  > sample_list.txt
for i in $(seq 1 22); do
  echo $i
  mkdir roh/chr$i
  cat sample_list.txt | xargs -I {} -P 16 sh -c "echo {}; bcftools roh --samples {} --threads 2 -G30 vcf/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz > roh/chr$i/{}_roh.txt; pigz -p2 roh/chr$i/{}_roh.txt"
done


# what we should do first is limit the VCFs to the sites that are present in the EGA haplotypecaller data... 
for i in $(seq 1 22); do
  echo "STARTING $i"
  bcftools annotate --rename-chrs rename_chrs.txt ~/data/p6/EGA_haplotypecaller/chr"$i"/chr"$i".vcf.gz -Oz  --threads 8 > EGA_chr"$i".vcf.gz
  tabix EGA_chr"$i".vcf.gz
  bcftools isec -p vcf_intersect_EGA/chr"$i" -Oz -c all -n=2 \
    vcf/ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
    EGA_chr"$i".vcf.gz
  rm EGA_chr"$i".vcf.gz
done

# this has the data we want vcf_intersect_EGA/chr1/0000.vcf.gz
# and I verified that it contains the right sites (1KG subsetted to the EGA whole exome)
for i in $(seq 1 22); do
  mv vcf_intersect_EGA/chr$i/0000.vcf.gz vcf_intersect_EGA/chr$i.vcf.gz
  mv vcf_intersect_EGA/chr$i/0000.vcf.gz.tbi vcf_intersect_EGA/chr$i.vcf.gz.tbi
done

# calculate ROH with the intersected data
for i in $(seq 8 22); do
  echo $i
  mkdir roh_intersect_EGA/chr$i
  cat sample_list.txt | xargs -I {} -P 16 sh -c "echo {}; bcftools roh --samples {} --threads 2 -G30 vcf_intersect_EGA/chr$i.vcf.gz > roh_intersect_EGA/chr$i/{}_roh.txt; pigz -p2 roh_intersect_EGA/chr$i/{}_roh.txt"
done

# now let's look at the patient data and see if we can to anything with the germline calls
# so probably want to run haplotypecaller on these, because before it was just mutect
# that's going to be a pita, because this data is completely unorganized. woot. 
outdir=/path/to/data/p6/CTP_201_haplotypecaller
cd $outdir
intervals=/path/to/data/references/AML_261_IDT/targets_w100_hg19.st.bed
ref=/path/to/data/references/human_g1k_v37_fixedChr.fasta
inbam=/home/user/box/shared_data/CTP_201/bam_all/DS-363463_012-201_Dx_PG_R1R2.bwa.mem.st.rmdup.bam
s="012-201_Dx_PG"
outgvcf=/path/to/data/p6/CTP_201_haplotypecaller/05_haplotypecaller/"$s".g.vcf.gz

###################################################################################################################################
############################# 2022-07-28 ##########################################################################################
###################################################################################################################################
# checking Y chromosome alignments with all EGA data 
cd /path/to/data/EGA_data/01_prepare_bam/recalibrate
ls *.bam | xargs -I {} -P 8 sh -c "echo {}; seqkit bam -C {} 2>&1 | cut -f 1,2 > /path/to/data/p6/EGA_Y_chromosome/seqkit_bam_out/{}.txt"

###################################################################################################################################
############################# 2022-07-30 ##########################################################################################
###################################################################################################################################
# trying alternative ROH detection tools
# AUDACITY is one that came up
cd /path/to/data/software/audacity
chr=9
vcf_zipped=/path/to/data/EGA_data_whole_exome/06_GDB/chr"$chr".vcf.gz
vcf=/path/to/data/EGA_data_whole_exome/06_GDB/chr"$chr".vcf
zcat $vcf_zipped > $vcf
outdir=/path/to/data/p6/EGA_haplotypecaller/audacity/chr"$chr"
AUDACITYPrepare.pl -I $vcf -F false -O $outdir -L chr"$chr"
AUDACITYAnalyze.pl -I "$outdir" -F "$outdir"/chr"$chr".HetCount.txt -O "$outdir" -R2 value -0.001 -AS hg19 -L chr"$chr"
# this software sucks and I'm never going to use it again. 

###################################################################################################################################
############################# 2022-08-07 ##########################################################################################
###################################################################################################################################
# getting back into some of this research. I need to run bcftools ROH 
# with the annotation and genetic map feature to ensure the results are legitimate. 
# should also look into filtering the results of the variant calling more, 
# I think the germline calls are not well filtered at present. 

cd /path/to/data/p6/EGA_haplotypecaller
for i in $(seq 1 22); do
  echo "#######################################################"
  echo "$i"
  echo "#######################################################"
  bcftools annotate --rename-chrs rename_chrs.txt chr"$i"/chr"$i".vcf.gz | \
    bcftools annotate -c CHROM,POS,REF,ALT,AF1KG \
    -h /path/to/data/p6/1000GP_data/1000GP-AFs/AFs.tab.gz.hdr \
    -a /path/to/data/p6/1000GP_data/1000GP-AFs/AFs_MOD.tab.gz -Oz --threads 16 > roh_annotated/chr"$i".vcf.gz

  bcftools roh --AF-tag AF1KG -M 100 \
    -m /path/to/data/p6/1000GP_data/genetic-map/genetic_map_chr"$i"_combined_b37.txt \
    -o roh_annotated/roh_chr"$i".txt \
    roh_annotated/chr"$i".vcf.gz
  grep RG roh_annotated/roh_chr"$i".txt > roh_annotated/roh_chr"$i"_RG.txt
  pigz roh_annotated/roh_chr"$i".txt
done


###################################################################################################################################
############################# 2022-08-30 ##########################################################################################
###################################################################################################################################
# have to re-call all bams 
# Repeat this for 102, 101 data as well
ls *.bam | xargs -I {} -P 8 sh -c "echo {}; bcftools mpileup -d 99999 -Ov -R ~/data/p6/bedfiles/p6_pattern_nochr.bed -f ~/data/references/human_g1k_v37_no_chr.fasta {} | bcftools call -c - > /path/to/data/p6/investigate_sample_data/CTP_201/call_bam/{}.vcf"
cd /path/to/data/p6/investigate_sample_data/CTP_201/call_bam/
ls *.vcf | xargs -I {} -P 16 sh -c "echo {}; bgzip {}; tabix {}.gz"
bcftools merge *.vcf.gz -Oz --threads 8 > merged.vcf.gz
plink2 --vcf merged.vcf.gz --make-bfile --out merged --max-alleles 2

#######################################################################
# Better to do it with GATK haplotypecaller 
# for each trial, parse out the reads around the P6 and GGCC variants from bams on box
cd /path/to/cloud_storage/shared_data/CTP_102/bam_all_20220819/
ls DS-38786*.bam | xargs -I {} -P 8 sh -c "echo {}; samtools view {} --regions-file ~/local_data/p6/p6_GGCC_targets_nochr.bed -O bam > /path/to/local_data/p6/investigate_sample_data/CTP_102/germline_call_bam/bam_subset/{}"
# then run the modified germline variant calling, which skips the bqsr and annotations
cd /path/to/local_data/p6/investigate_sample_data
snakemake -s ~/projects/project/call_bam_GATK/call_bam_GATK_germline.snakefile --configfile config_call_bam_GATK_germline_CTP_102.yaml --jobs 24 --cores 24 --rerun-incomplete 
# convert the variant files with plink
plink2 --vcf CTP_102/germline_call_bam/06_GDB/chr9/chr9.vcf --make-bfile --max-alleles 2 --out CTP_102/germline_call_bam/06_GDB/chr9/chr9

