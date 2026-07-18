## AWS EC2 setup notes (historical)
# These are rough notes from the original AWS/S3-based execution model, kept for reference
# only. Credentials have been stripped. AWS support will be rewritten from scratch as a
# separate effort - treat this as background, not a working script.

ssh -i ~/aws/keys/<key-name>.pem ubuntu@<ec2-host>

sudo apt update
sudo apt install awscli -y


sudo mkfs -t ext4 /dev/nvme3n1
sudo mkdir /fasttmp
sudo mount /dev/nvme3n1 /fasttmp
sudo chmod -R 777 /fasttmp

sudo mkfs -t ext4 /dev/nvme1n1
sudo mkdir ~/references
sudo mount /dev/nvme1n1 ~/references
sudo chmod -R 777 ~/references
mkdir ~/references/data-bucket
ln -s /home/ubuntu/references/data-bucket /home/ubuntu/

sudo mkfs -t ext4 /dev/nvme4n1
sudo mkdir ~/output
sudo mount /dev/nvme4n1 ~/output
sudo chmod -R 777 ~/output

# maybe do something about setting TMPDIR to /fasttmp on startup?
echo "export TMPDIR=/fasttmp" >> ~/.bashrc
echo "export TMP=/fasttmp" >> ~/.bashrc
echo "export TEMP=/fasttmp" >> ~/.bashrc
source ~/.bashrc
echo "  " > ~/.tmux.conf


sudo mkfs -t ext4 /dev/nvme2n1
sudo mkdir /swap2
sudo mount /dev/nvme2n1 /swap2

sudo fallocate -l 250G ~/swap
sudo chmod 600 ~/swap
sudo mkswap ~/swap
sudo swapon ~/swap

aws s3 sync s3://<your-bucket>/references references &
aws s3 sync s3://<your-bucket>/data data &
aws s3 sync s3://<your-bucket>/projects projects &

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh -b
./miniconda3/condabin/conda init
source ~/.bashrc
conda install mamba -n base -c conda-forge -y
mamba install -c conda-forge -c bioconda gatk4 snakemake=6 bwa samtools bedtools bcftools pyega3 picard -y


snakemake -s ~/projects/pipeline/call_bam_gatk/AWS_call_bam_gatk.snakefile --configfile config_call_bam_GATK.yaml --jobs 64 --cores 64 -k

# [REDACTED] EGA/AWS credentials used to live here in plaintext - removed.
# Store credentials in ~/.aws/credentials and a separate EGA credential_file.json instead,
# never inline in a script.

time pyega3 -c 30 fetch -M -1 -W 2 --format BAM EGAF00000064163

cat not_downloaded.txt | xargs -P 16 -I {} bash dl.sh {}
# test if the bam exists on AWS yet
# this is the dl.sh script referenced above
if aws s3 ls s3://<your-bucket>/data/"$1"/ | grep bam; then
        echo "$1 EXISTS ALREADY! NOT DOWNLOADING"
else
        echo "$1 STARTING "
        pyega3 -c 60 fetch -M -1 -W 2 --format BAM $1 && aws s3 cp --recursive "$1" s3://<your-bucket>/data/"$1" && rm -r "$1"
        while [ $? -ne 0 ]; do
            echo "RESTARTING $1"
            pyega3 -c 60 fetch -M -1 -W 2 --format BAM $1 && aws s3 cp --recursive "$1" s3://<your-bucket>/data/"$1" && rm -r "$1"
            sleep 1
        done
fi
