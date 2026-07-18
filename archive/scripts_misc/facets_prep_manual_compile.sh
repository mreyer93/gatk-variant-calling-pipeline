# compile snp-pileup based on where htslib is installed in miniconda
g++ -std=c++11 -I/path/to/data/miniconda3/include/htslib/ \
    snp-pileup.cpp  -L/path/to/data/miniconda3/lib \
    -lhts -Wl,-rpath=/path/to/data/miniconda3/lib \
    -o snp-pileup                                 


snp-pileup <vcf file> <output file> <sequence files...>

# try this on some CTP-102 samples 
# VCF file is dbsnp from NCBI
pt=001-101
vcff=/path/to/data/references/dbsnp/human_9606_b151_GRCh37p13/00-common_all.vcf.gz
outf=006-006.csv.gz
infs="/path/to/data/ctp_101_calling/call_bam_new/output/01_prepare_bam/recalibrate/006-006_Dx_HS.fixchr.bam /path/to/data/ctp_101_calling/call_bam_new/output/01_prepare_bam/recalibrate/006-006_Dx_BM.fixchr.bam"
# normal-tumor order!
rm -f $outf.gz && snp-pileup  -g -q15 -Q20 -P100 -r25,0 -v -d 10000 $vcff $outf $infs

# need to run this for all normal-tumor samples
