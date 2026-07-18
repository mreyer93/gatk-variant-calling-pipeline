zcat 1000G_omni2.5.b37.vcf.gz | \
sed -e 's/contig=<ID=1,/contig=<ID=chr1,/' | sed -e 's/contig=<ID=2,/contig=<ID=chr2,/' | \
sed -e 's/contig=<ID=3,/contig=<ID=chr3,/' | sed -e 's/contig=<ID=4,/contig=<ID=chr4,/' | \
sed -e 's/contig=<ID=5,/contig=<ID=chr5,/' | sed -e 's/contig=<ID=6,/contig=<ID=chr6,/' | \
sed -e 's/contig=<ID=7,/contig=<ID=chr7,/' | sed -e 's/contig=<ID=8,/contig=<ID=chr8,/' | \
sed -e 's/contig=<ID=9,/contig=<ID=chr9,/' | sed -e 's/contig=<ID=10,/contig=<ID=chr10,/' | \
sed -e 's/contig=<ID=11,/contig=<ID=chr11,/' | sed -e 's/contig=<ID=12,/contig=<ID=chr12,/' | \
sed -e 's/contig=<ID=13,/contig=<ID=chr13,/' | sed -e 's/contig=<ID=14,/contig=<ID=chr14,/' | \
sed -e 's/contig=<ID=15,/contig=<ID=chr15,/' | sed -e 's/contig=<ID=16,/contig=<ID=chr16,/' | \
sed -e 's/contig=<ID=17,/contig=<ID=chr17,/' | sed -e 's/contig=<ID=18,/contig=<ID=chr18,/' | \
sed -e 's/contig=<ID=19,/contig=<ID=chr19,/' | sed -e 's/contig=<ID=20,/contig=<ID=chr20,/' | \
sed -e 's/contig=<ID=21,/contig=<ID=chr21,/' | sed -e 's/contig=<ID=22,/contig=<ID=chr22,/' | \
sed -e 's/contig=<ID=X,/contig=<ID=chrX,/' | sed -e 's/contig=<ID=Y,/contig=<ID=chrY,/' | \
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
bgzip  > 1000G_omni2.5.b37.chr.vcf.gz
gatk IndexFeatureFile _I 1000G_omni2.5.b37.chr.vcf.gz


zcat 1000G_phase1.snps.high_confidence.b37.vcf.gz | \
sed -e 's/contig=<ID=1,/contig=<ID=chr1,/' | sed -e 's/contig=<ID=2,/contig=<ID=chr2,/' | \
sed -e 's/contig=<ID=3,/contig=<ID=chr3,/' | sed -e 's/contig=<ID=4,/contig=<ID=chr4,/' | \
sed -e 's/contig=<ID=5,/contig=<ID=chr5,/' | sed -e 's/contig=<ID=6,/contig=<ID=chr6,/' | \
sed -e 's/contig=<ID=7,/contig=<ID=chr7,/' | sed -e 's/contig=<ID=8,/contig=<ID=chr8,/' | \
sed -e 's/contig=<ID=9,/contig=<ID=chr9,/' | sed -e 's/contig=<ID=10,/contig=<ID=chr10,/' | \
sed -e 's/contig=<ID=11,/contig=<ID=chr11,/' | sed -e 's/contig=<ID=12,/contig=<ID=chr12,/' | \
sed -e 's/contig=<ID=13,/contig=<ID=chr13,/' | sed -e 's/contig=<ID=14,/contig=<ID=chr14,/' | \
sed -e 's/contig=<ID=15,/contig=<ID=chr15,/' | sed -e 's/contig=<ID=16,/contig=<ID=chr16,/' | \
sed -e 's/contig=<ID=17,/contig=<ID=chr17,/' | sed -e 's/contig=<ID=18,/contig=<ID=chr18,/' | \
sed -e 's/contig=<ID=19,/contig=<ID=chr19,/' | sed -e 's/contig=<ID=20,/contig=<ID=chr20,/' | \
sed -e 's/contig=<ID=21,/contig=<ID=chr21,/' | sed -e 's/contig=<ID=22,/contig=<ID=chr22,/' | \
sed -e 's/contig=<ID=X,/contig=<ID=chrX,/' | sed -e 's/contig=<ID=Y,/contig=<ID=chrY,/' | \
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
bgzip  > 1000G_phase1.snps.high_confidence.b37.chr.vcf.gz
gatk IndexFeatureFile -I 1000G_phase1.snps.high_confidence.b37.chr.vcf.gz


zcat dbsnp_138.b37.vcf.gz | \
sed -e 's/contig=<ID=1,/contig=<ID=chr1,/' | sed -e 's/contig=<ID=2,/contig=<ID=chr2,/' | \
sed -e 's/contig=<ID=3,/contig=<ID=chr3,/' | sed -e 's/contig=<ID=4,/contig=<ID=chr4,/' | \
sed -e 's/contig=<ID=5,/contig=<ID=chr5,/' | sed -e 's/contig=<ID=6,/contig=<ID=chr6,/' | \
sed -e 's/contig=<ID=7,/contig=<ID=chr7,/' | sed -e 's/contig=<ID=8,/contig=<ID=chr8,/' | \
sed -e 's/contig=<ID=9,/contig=<ID=chr9,/' | sed -e 's/contig=<ID=10,/contig=<ID=chr10,/' | \
sed -e 's/contig=<ID=11,/contig=<ID=chr11,/' | sed -e 's/contig=<ID=12,/contig=<ID=chr12,/' | \
sed -e 's/contig=<ID=13,/contig=<ID=chr13,/' | sed -e 's/contig=<ID=14,/contig=<ID=chr14,/' | \
sed -e 's/contig=<ID=15,/contig=<ID=chr15,/' | sed -e 's/contig=<ID=16,/contig=<ID=chr16,/' | \
sed -e 's/contig=<ID=17,/contig=<ID=chr17,/' | sed -e 's/contig=<ID=18,/contig=<ID=chr18,/' | \
sed -e 's/contig=<ID=19,/contig=<ID=chr19,/' | sed -e 's/contig=<ID=20,/contig=<ID=chr20,/' | \
sed -e 's/contig=<ID=21,/contig=<ID=chr21,/' | sed -e 's/contig=<ID=22,/contig=<ID=chr22,/' | \
sed -e 's/contig=<ID=X,/contig=<ID=chrX,/' | sed -e 's/contig=<ID=Y,/contig=<ID=chrY,/' | \
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
bgzip  > dbsnp_138.b37.chr.vcf.gz
gatk IndexFeatureFile -I dbsnp_138.b37.chr.vcf.gz

zcat hapmap_3.3.b37.vcf.gz | \
sed -e 's/contig=<ID=1,/contig=<ID=chr1,/' | sed -e 's/contig=<ID=2,/contig=<ID=chr2,/' | \
sed -e 's/contig=<ID=3,/contig=<ID=chr3,/' | sed -e 's/contig=<ID=4,/contig=<ID=chr4,/' | \
sed -e 's/contig=<ID=5,/contig=<ID=chr5,/' | sed -e 's/contig=<ID=6,/contig=<ID=chr6,/' | \
sed -e 's/contig=<ID=7,/contig=<ID=chr7,/' | sed -e 's/contig=<ID=8,/contig=<ID=chr8,/' | \
sed -e 's/contig=<ID=9,/contig=<ID=chr9,/' | sed -e 's/contig=<ID=10,/contig=<ID=chr10,/' | \
sed -e 's/contig=<ID=11,/contig=<ID=chr11,/' | sed -e 's/contig=<ID=12,/contig=<ID=chr12,/' | \
sed -e 's/contig=<ID=13,/contig=<ID=chr13,/' | sed -e 's/contig=<ID=14,/contig=<ID=chr14,/' | \
sed -e 's/contig=<ID=15,/contig=<ID=chr15,/' | sed -e 's/contig=<ID=16,/contig=<ID=chr16,/' | \
sed -e 's/contig=<ID=17,/contig=<ID=chr17,/' | sed -e 's/contig=<ID=18,/contig=<ID=chr18,/' | \
sed -e 's/contig=<ID=19,/contig=<ID=chr19,/' | sed -e 's/contig=<ID=20,/contig=<ID=chr20,/' | \
sed -e 's/contig=<ID=21,/contig=<ID=chr21,/' | sed -e 's/contig=<ID=22,/contig=<ID=chr22,/' | \
sed -e 's/contig=<ID=X,/contig=<ID=chrX,/' | sed -e 's/contig=<ID=Y,/contig=<ID=chrY,/' | \
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
bgzip  > hapmap_3.3.b37.chr.vcf.gz
gatk IndexFeatureFile -I hapmap_3.3.b37.chr.vcf.gz