# Variant calling in sample-vs-reference mode is done for whole exome
# Variant calling in tumor-vs-normal mode NOT DONE for whole exome
# Variant calling germline NOT DONE in whole exome

cd /home/ubuntu/ega_recall/
export SENTIEON_LICENSE=/home/ubuntu/Formic_Labs_eval.lic
SENTIEON_INSTALL_DIR=/home/ubuntu/sentieon-genomics-202112.06
NT=16

# Update with the location of the reference data files
FASTA_DIR="/home/ubuntu/references/"
FASTA="$FASTA_DIR/human_g1k_v37_fixedChr.fasta"
KNOWN_DBSNP="$FASTA_DIR/dbsnp/human_9606_b151_GRCh37p13/00-common_all.vcf.gz"
INTERVAL_FILE="$FASTA_DIR/sure_select_V5/S04380110_Padded.bed"
BAM_DIR=/home/ubuntu/ega_recall/bam2/
LOG_DIR=/home/ubuntu/ega_recall/log/
OUTDIR_GVCF=/home/ubuntu/ega_recall/sentieon_GVCF_output
OUTDIR_HAPLOTYPER=/home/ubuntu/ega_recall/sentieon_Haplotyper_output


# !!! PD6606b ERROR!!!

for SAMPLE in $(cat /home/ubuntu/ega_recall/sample_list.txt); do
    echo "$SAMPLE"
    BAM=$BAM_DIR/$SAMPLE.fixchr.bam

    if [ -f $OUTDIR_HAPLOTYPER/$SAMPLE.vcf.gz ] && [ -f $OUTDIR_GVCF/$SAMPLE.g.vcf.gz ]; then
        continue
    fi
    rm -f $LOG_DIR/stage_$SAMPLE.log
    rm -f $LOG_DIR/Haplotyper_$SAMPLE.log
    rm -f $LOG_DIR/gvcf_$SAMPLE.log

    echo "STAGING $SAMPLE"
    (time rclone copy box:shared_data/EGA_data/01_prepare_bam/recalibrate/$SAMPLE.fixchr.bam $BAM_DIR) 2>&1 | tee -a $LOG_DIR/stage_$SAMPLE.log
    (time rclone copy box:shared_data/EGA_data/01_prepare_bam/recalibrate/$SAMPLE.fixchr.bam.bai $BAM_DIR) 2>&1 | tee -a $LOG_DIR/stage_$SAMPLE.log

    echo "HAPLOTYPER $SAMPLE"
    $SENTIEON_INSTALL_DIR/bin/sentieon driver \
        -r $FASTA \
        ${INTERVAL_FILE:+--interval $INTERVAL_FILE} \
        -t $NT \
        -i $BAM \
        --algo Haplotyper \
        -d $KNOWN_DBSNP \
        --emit_conf=30 \
        --call_conf=30 \
        $OUTDIR_HAPLOTYPER/$SAMPLE.vcf.gz 2>&1 | tee -a $LOG_DIR/Haplotyper_$SAMPLE.log

    echo "GVCF $SAMPLE"
    $SENTIEON_INSTALL_DIR/bin/sentieon driver \
        -r $FASTA \
        ${INTERVAL_FILE:+--interval $INTERVAL_FILE} \
        -t $NT \
        -i $BAM \
        --algo Haplotyper \
        -d $KNOWN_DBSNP \
        --emit_mode gvcf \
        $OUTDIR_GVCF/$SAMPLE.g.vcf.gz 2>&1 | tee -a $LOG_DIR/gvcf_$SAMPLE.log

    rm -f $BAM_DIR/$SAMPLE*
done

# --emit_mode gvcf
# gvcf: emits additional information required for joint calling. This option is required if you want to perform joint calling using the GVCFtyper algorithm

$SENTIEON_INSTALL_DIR/bin/sentieon driver \
    -r $FASTA \
    --algo GVCFtyper \
    -d $KNOWN_DBSNP \
    --emit_conf=30 \
    --call_conf=30 \
    output-joint.vcf.gz \
    $OUTDIR_GVCF/*.g.vcf.gz

# Norm and LD prune
bcftools norm \
    --multiallelics -both \
    --output-type z \
     --threads 15 \
    --output output-joint_norm.vcf.gz \
    output-joint.vcf.gz

bcftools +prune \
    -m 0.2 \
    --output-type z \
    --output output-joint_norm_LD.vcf.gz \
    output-joint_norm.vcf.gz

# convert to tsv
bcftools query \
    --print-header \
    -f "%CHROM\t%POS\t%TYPE\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/DP\t%INFO/AF\t%INFO/AN\t%INFO/AC[\t%GT\t%AD\t%DP\t%GQ\t%PL]\n" \
     output-joint_norm_LD.vcf.gz \
    -o output-joint_norm_LD.tsv


BAM_DIR=/home/ubuntu/ssd1/bam/
tail -n +310 /home/ubuntu/ega_recall/sample_list.txt | head -n 10 | xargs -P 16 -I {} sh -c "
    echo {}
    echo STAGING {}
    rclone copy -P box:shared_data/EGA_data/01_prepare_bam/recalibrate/{}.fixchr.bam $BAM_DIR
    rclone copy -P box:shared_data/EGA_data/01_prepare_bam/recalibrate/{}.fixchr.bam.bai $BAM_DIR
    "

BAM_DIR=/home/ubuntu/ssd2/bam/
tail -n 18 /home/ubuntu/ega_recall/sample_list.txt  | xargs -P 8 -I {} sh -c "
    echo {}
    echo STAGING {}
    rclone copy -P box:shared_data/EGA_data/01_prepare_bam/recalibrate/{}.fixchr.bam $BAM_DIR
    rclone copy -P box:shared_data/EGA_data/01_prepare_bam/recalibrate/{}.fixchr.bam.bai $BAM_DIR
    "

# Doing annotation with funcotator, locally
cd /path/to/local_data/EGA_recall/output-joint
gatk Funcotator \
    --variant output-joint_norm.vcf \
    --reference /path/to/local_data/references/human_g1k_v37_fixedChr.fasta \
    --ref-version hg19 \
    --data-sources-path /path/to/local_data/references/funcotator_dataSources.v1.7.20200521g \
    --output output-joint_norm_func.vcf \
    --output-file-format VCF