cd genomes; wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz & wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

split-pipe \
	--mode mkref \
	--genome_name hg38 \
	--fasta ~/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
	--genes ~/genomes/Homo_sapiens.GRCh38.109.gtf.gz \
	--output_dir ~/genomes/hg38

MV411_1_SublibraryIndex1_S2_R1_001.fastq.gz
MV411_1_SublibraryIndex1_S2_R2_001.fastq.gz
MV411_2_SublibraryIndex2_S3_R1_001.fastq.gz
MV411_2_SublibraryIndex2_S3_R2_001.fastq.gz

split-pipe \
    --mode all \
    --chemistry v2 \
    --genome_dir ~/genomes/hg38 \
    --output_dir parse_analysis/sl1 \
    --nthreads 16 \
    --fq1 parse_data/MV411_1_SublibraryIndex1_S2_R1_001.fastq.gz

split-pipe \
    --mode all \
    --chemistry v2 \
    --genome_dir ~/genomes/hg38 \
    --output_dir parse_analysis/sl2 \
    --nthreads 16 \
    --fq1 parse_data/MV411_2_SublibraryIndex2_S3_R1_001.fastq.gz

# try by sample: ACTUALLY WANT TO DO IT THIS WAY SO WE CAN ANALYZE EACH WELL SEPARATELY
split-pipe \
    --mode all \
    --chemistry v2 \
    --genome_dir ~/genomes/hg38 \
    --output_dir ~/parse_analysis/sl1_bysample \
    --nthreads 16 \
    --fq1 ~/parse_data/MV411_1_SublibraryIndex1_S2_R1_001.fastq.gz \
    --sample A1 A1 \
    --sample A2 A2 \
    --sample A3 A3 \
    --sample A4 A4 \
    --sample A5 A5 \
    --sample A6 A6 \
    --sample A7 A7 \
    --sample A8 A8 \
    --sample A9 A9 \
    --sample A10 A10 \
    --sample A11 A11 \
    --sample A12 A12
split-pipe \
    --mode all \
    --chemistry v2 \
    --genome_dir ~/genomes/hg38 \
    --output_dir ~/parse_analysis/sl2_bysample \
    --nthreads 16 \
    --fq1 ~/parse_data/MV411_2_SublibraryIndex2_S3_R1_001.fastq.gz \
    --sample A1 A1 \
    --sample A2 A2 \
    --sample A3 A3 \
    --sample A4 A4 \
    --sample A5 A5 \
    --sample A6 A6 \
    --sample A7 A7 \
    --sample A8 A8 \
    --sample A9 A9 \
    --sample A10 A10 \
    --sample A11 A11 \
    --sample A12 A12


# combine sublibraries
split-pipe --mode comb \
	--sublibraries ~/parse_analysis/sl1_bysample ~/parse_analysis/sl2_bysample \
	--output_dir ~/parse_analysis/combined 

# get to h5ad
cd ~/parse_analysis
echo "ana_save_anndata True" > parfile.txt

split-pipe --mode ana \
    --nthreads 16 \
    --genome_dir ~/genomes/hg38 \
    --output_dir sl1 \
    --parfile parfile.txt \
    --chemistry v2 \
    --kit WT_mini \
    --sample A1 A1 \
    --sample A2 A2 \
    --sample A3 A3 \
    --sample A4 A4 \
    --sample A5 A5 \
    --sample A6 A6 \
    --sample A7 A7 \
    --sample A8 A8 \
    --sample A9 A9 \
    --sample A10 A10 \
    --sample A11 A11 \
    --sample A12 A12

# process each well separately
split-pipe --mode ana \
    --nthreads 16 \
    --genome_dir ~/genomes/hg38 \
    --output_dir combined \
    --parfile parfile.txt \
    --chemistry v2 \
    --kit WT_mini \
    --sample A1 A1 \
    --sample A2 A2 \
    --sample A3 A3 \
    --sample A4 A4 \
    --sample A5 A5 \
    --sample A6 A6 \
    --sample A7 A7 \
    --sample A8 A8 \
    --sample A9 A9 \
    --sample A10 A10 \
    --sample A11 A11 \
    --sample A12 A12