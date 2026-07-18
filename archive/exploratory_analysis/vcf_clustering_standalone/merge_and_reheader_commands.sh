# commands to convert and merge vcf
# first bgzip and index with bcftools index
# merge them all into one

# bcftools merge -f PASS  -O z  -o  ../merged_filtered_no0.vcf.gz --threads 32 *filtered.vcf.gz
# bcftools merge -f PASS -0 -O z  -o  ../merged_filtered.vcf.gz --threads 32 *filtered.vcf.gz

# need to remove the AS_SB_TABLE field because of problems with seqvcf2gds
# AS_SB_TABLE=60,2|48,2;


bcftools view -h merged_filtered.vcf.gz | grep -v "AS_SB_TABLE" > merged_filtered_fixed.header
zcat merged_filtered.vcf.gz | sed "s|AS_SB_TABLE=[^ECNT]*||g"  | bgzip | \
    bcftools reheader -h merged_filtered_fixed.header  -o merged_filtered_fixed.vcf.gz --threads 32 -
mv merged_filtered_fixed.vcf.gz merged_filtered.vcf.gz
rm merged_filtered_fixed.header

bcftools view -h merged_filtered_no0.vcf.gz | grep -v "AS_SB_TABLE" > merged_filtered_no0_fixed.header
zcat merged_filtered_no0.vcf.gz | sed "s|AS_SB_TABLE=[^ECNT]*||g"  | bgzip | \
    bcftools reheader -h merged_filtered_no0_fixed.header  -o merged_filtered_no0_fixed.vcf.gz --threads 32 -
mv merged_filtered_no0_fixed.vcf.gz merged_filtered_no0.vcf.gz
rm merged_filtered_no0_fixed.header



# count number of reads in each aml targets coverage file
cut -f 6  PD6637a2.w100.cov | awk '{s+=$1} END {print s}'

rm -f total_reads_aml_targets.txt
for i in *.cov; do 
    echo $i
    echo -e "$i""\t"$(cut -f 6 $i | awk '{s+=$1} END {print s}')>> total_reads_aml_targets.txt
done
