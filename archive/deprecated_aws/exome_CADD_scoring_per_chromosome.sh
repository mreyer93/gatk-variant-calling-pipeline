# take unique variants from all samples
cd /path/to/data/EGA_data/02_variants_reference_exome
mkdir unique_exome
zcat filtered/*.gz | cut -f 1-5 | sort | uniq > unique_exome/variants.vcf
cd unique_exome

# separate out into chromosomes
grep -P "^chr1\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_01.sorted.vcf.gz
grep -P "^chr2\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_02.sorted.vcf.gz
grep -P "^chr3\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_03.sorted.vcf.gz
grep -P "^chr4\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_04.sorted.vcf.gz
grep -P "^chr5\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_05.sorted.vcf.gz
grep -P "^chr6\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_06.sorted.vcf.gz
grep -P "^chr7\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_07.sorted.vcf.gz
grep -P "^chr8\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_08.sorted.vcf.gz
grep -P "^chr9\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_09.sorted.vcf.gz
grep -P "^chr10\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_10.sorted.vcf.gz
grep -P "^chr11\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_11.sorted.vcf.gz
grep -P "^chr12\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_12.sorted.vcf.gz
grep -P "^chr13\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_13.sorted.vcf.gz
grep -P "^chr14\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_14.sorted.vcf.gz
grep -P "^chr15\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_15.sorted.vcf.gz
grep -P "^chr16\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_16.sorted.vcf.gz
grep -P "^chr17\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_17.sorted.vcf.gz
grep -P "^chr18\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_18.sorted.vcf.gz
grep -P "^chr19\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_19.sorted.vcf.gz
grep -P "^chr20\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_20.sorted.vcf.gz
grep -P "^chr21\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_21.sorted.vcf.gz
grep -P "^chr22\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_22.sorted.vcf.gz
grep -P "^chrX\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_X.sorted.vcf.gz
grep -P "^chrY\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_Y.sorted.vcf.gz
grep -P "^chrM\t" variants.vcf | sed -s "s/^chr//g" | sort -k 2,2 -n | bgzip > variants_M.sorted.vcf.gz

# Split up based on file line numbers and compute that way
# parallelize across many sets
cd /path/to/data/EGA_data/02_variants_reference_exome/unique_exome
# chrs=( 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y )
chrs=( 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 )
n_lines=1000
n_processes=64
for chr in "${chrs[@]}"; do
    echo "STARTING " $chr
    outdir=variants_"$chr"_split
    mkdir -p "$outdir"
    cd "$outdir"
    inf=../variants_"$chr".sorted.vcf.gz
    zcat "$inf" | split -l $n_lines --additional-suffix ".vcf"
    ls x*.vcf | xargs -P $n_processes -I {} sh -c "~/cadd_xargs.sh {}"
    zcat x*.tsv.gz | grep -v "#" > cadd_combined.tsv
    paste <(cut -f 1-4 cadd_combined.tsv | tr "\t" "_") cadd_combined.tsv | sort -k1,1 > cadd_combined_join.tsv
    rm cadd_combined.tsv
    cd /path/to/data/EGA_data/02_variants_reference_exome/unique_exome
done

# combine all of them - have to sort anyway so not going to follow chr order
cd /path/to/data/EGA_data/02_variants_reference_exome/unique_exome
cat variants_*_split/cadd_combined_join.tsv | sort -k 1b,1 > allchr_combined.txt

# do a lookup / join with the actual data
for f in /path/to/data/EGA_data/02_variants_reference_exome/filtered/*.vcf.gz; do
zcat 1711STDY5270146_filtered.vcf.gz | grep -v "#" | sed "s/^chr//g" | cut -f 1,2,4,5 | tr "\t" "_"  | sort -k 1,1 | join - allchr_combined.txt > 

cut -f 1,2,4,5 1711STDY5270146.tmp.vcf  | tr "\t" "_" | sort -k 1,1 > 1711STDY5270146.lookup
join -t "\t" 1711STDY5270146.lookup allchr_combined.txt > 1711STDY5270146.cadd.tsv


# do across everything 
ls /path/to/data/EGA_data/02_variants_reference_exome/filtered/*.vcf.gz | cut -f 8 -d "/" | sed "s/.vcf.gz//g" > sample_lookup_list.txt                         
cat sample_lookup_list.txt | xargs -I {} -P 60 sh -c 'echo {}; zcat /path/to/data/EGA_data/02_variants_reference_exome/filtered/{}.vcf.gz | grep -v "#" | sed "s/^chr//g" | cut -f 1,2,4,5 | tr "\t" "_"  | sort -k 1b,1 | join - allchr_combined.txt | tr " " "\t" > cadd_lookup_results/{}.cadd.tsv'


# then this needs to feed back in the snakemake part of the script, where it gets separated into the simple script
outf=/path/to/data/EGA_data/02_variants_reference_exome/annotate_CADD/simple_table/1711STDY5270146_filtered_CADD_simple.vcf

cut -f 1,2,3,4,5,6,10,21,116 /path/to/data/EGA_data/02_variants_reference_exome/annotate_CADD/header.tsv > $outf
cut -f 2,3,4,5,6,7,11,22,117 cadd_lookup_results/1711STDY5270146_filtered.cadd.tsv | sort -V -k1,1 -k2,2  >> $outf
gzip $outf

# across ALL
cd /path/to/data/EGA_data/02_variants_reference_exome/unique_exome
cat sample_lookup_list.txt | xargs -I {} -P 60 sh -c 'echo {}; cut -f 1,2,3,4,5,6,10,21,116 /path/to/data/EGA_data/02_variants_reference_exome/annotate_CADD/header.tsv > /path/to/data/EGA_data/02_variants_reference_exome/annotate_CADD/simple_table/{}_CADD_simple.vcf; cut -f 2,3,4,5,6,7,11,22,117 cadd_lookup_results/{}.cadd.tsv | sort -V -k1,1 -k2,2  >> /path/to/data/EGA_data/02_variants_reference_exome/annotate_CADD/simple_table/{}_CADD_simple.vcf; gzip -f /path/to/data/EGA_data/02_variants_reference_exome/annotate_CADD/simple_table/{}_CADD_simple.vcf'