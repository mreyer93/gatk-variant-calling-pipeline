# transforming gnomAD csv file into a VCF that can be used with CADD
# only want a few simple columns
# CHROM	POS	ID	REF	ALT

inf="gnomAD_BOD1L1.csv"
outf="gnomAD_BOD1L1.vcf"
outf_cadd="gnomAD_BOD1L1_cadd_orig.vcf.gz"
outf_cadd2="gnomAD_BOD1L1_cadd.vcf"
outf_cadd3="gnomAD_BOD1L1_cadd_subset.vcf"
inf="gnomAD_MAP1B.csv"
outf="gnomAD_MAP1B.vcf"
outf_cadd="gnomAD_MAP1B_cadd_orig.vcf.gz"
outf_cadd2="gnomAD_MAP1B_cadd.vcf"
outf_cadd3="gnomAD_MAP1B_cadd_subset.vcf"
echo -e "#CHROM\tPOS\tID\tREF\tALT" > $outf
tail -n +2 $inf | awk 'BEGIN {FS=",";OFS="\t"} {print $1,$2,$3,$4,$5}' >> $outf

CADD.sh -a -c 4 -g GRCh37 -o $outf_cadd $outf
zcat $outf_cadd | cut -f 1,2,3,4,5,6,10,21,116 > $outf_cadd2

# NEED TO SUBSET THIS TO VARIANTS IN THE AML TARGET SET
bgzip $outf_cadd2
tabix -p vcf $outf_cadd2.gz
tabix -R targets_w100_hg19.st.nochr.bed $outf_cadd2.gz > $outf_cadd3

# now process it in R
"Chrom", "Pos", "Ref", "Alt", "Type", "Length", "ConsDetail", "GeneName", "PHRED"
