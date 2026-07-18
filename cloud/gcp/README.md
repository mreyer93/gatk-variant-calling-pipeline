# Running the full pipeline on GCP

This is the "full" version of the pipeline (with CADD, no annotation shortcuts) for when
your local machine doesn't have the 400G+ of disk space CADD needs. It's one big Compute
Engine VM you start, use, and stop - not a production/autoscaling setup. That matches
occasional personal-research/consulting use better than a constantly-running cluster
would, and it's what these scripts automate.

Why GCP over AWS: reference-data hosting is a wash (both [Broad's GRCh38 bundle](https://registry.opendata.aws/broad-references/)
and the GATK best-practices resources are freely mirrored on both), but GATK/Broad's own
tooling and tutorials are traditionally GCP-first, and spot VMs are a simple, well-known
way to cut compute cost for exactly this kind of occasional batch workload.

## Prerequisites
- A GCP project with billing enabled and the Compute Engine API turned on.
- [`gcloud` CLI](https://cloud.google.com/sdk/docs/install) installed and authenticated (`gcloud init`, `gcloud auth login`).
- Edit `cloud/gcp/config.sh` and set `PROJECT_ID` (and adjust `ZONE`/`MACHINE_TYPE`/`DISK_SIZE_GB` if you want).

## One-time setup
```bash
./cloud/gcp/create_disk.sh   # creates the persistent data disk (references, CADD, outputs)
./cloud/gcp/start_vm.sh      # creates the VM, attaches the disk, formats/mounts it at /mnt/data
gcloud compute ssh gatk-pipeline-vm --project=<your-project> --zone=us-central1-a
```
Once logged into the VM, set up the pipeline the same way as the [local instructions](../../manual/requirements.md), just pointed at `/mnt/data`:
```bash
# clone the repo - do this interactively so your GitHub credentials never touch VM
# metadata/scripts. Install gh and `gh auth login`, or use an SSH deploy key.
git clone https://github.com/<you>/gatk-variant-calling-pipeline.git
cd gatk-variant-calling-pipeline

# conda/mamba
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p /mnt/data/miniforge3
source /mnt/data/miniforge3/etc/profile.d/conda.sh
mamba env create -f envs/environment.yml -n gatk-pipeline
conda activate gatk-pipeline

# references (full set, onto the data disk - this is the ~17G core set; add
# --with-alphamissense too if you also want it for germline)
./scripts/download_references.sh /mnt/data/references
gatk FuncotatorDataSourceDownloader --germline --hg38 --validate-integrity --extract-after-download
gatk FuncotatorDataSourceDownloader --somatic  --hg38 --validate-integrity --extract-after-download
cadd-install.sh   # answer "no" to installing a separate CADD conda env - you already have one
```
This step takes a while (400G+ for CADD alone) but only happens once - it all lives on
the persistent disk, not the VM, so it survives stopping/deleting the VM.

## Routine workflow
```bash
./cloud/gcp/start_vm.sh
gcloud compute ssh gatk-pipeline-vm --project=<your-project> --zone=us-central1-a
# on the VM:
source /mnt/data/miniforge3/etc/profile.d/conda.sh && conda activate gatk-pipeline
cd gatk-variant-calling-pipeline
snakemake -s call_bam_GATK/call_bam_GATK.snakefile --configfile /mnt/data/my_run/config_call_bam_GATK.yaml --use-conda --jobs 16 --cores 16 -k
# copy results off before stopping, e.g. to a GCS bucket:
gsutil -m cp -r /mnt/data/my_run/output gs://<your-bucket>/my_run/
# back on your laptop:
./cloud/gcp/stop_vm.sh
```
`--jobs`/`--cores` should match `MACHINE_TYPE` in `config.sh` (16 for the n2-standard-16 default).

### If the spot VM gets preempted mid-run
Spot VMs can be reclaimed by GCP at any time (that's why they're ~70% cheaper). If that
happens mid-pipeline, just `./cloud/gcp/start_vm.sh` again and re-run the exact same
`snakemake` command - Snakemake sees the outputs that already completed and picks up
where it left off, it doesn't restart from scratch.

## Cost
Rough estimate, `us-central1`, current as of when this was written (verify with the
[pricing calculator](https://cloud.google.com/products/calculator) - prices drift):
- **Compute**: `n2-standard-16` spot ≈ $0.23/hr. Only billed while the VM is running.
- **Storage**: 1000G `pd-balanced` ≈ $100/month, billed **whether the VM is running or
  not** - this is the dominant cost for occasional use, not compute.

For truly occasional use (a few runs a year), storage cost while idle can matter more
than compute. Two ways to handle that if $100/month sitting idle bothers you:
1. **Full teardown between uses** (below) and re-download references/re-run
  `cadd-install.sh` each time - saves the storage cost, costs you the download time
  (hours) on your next run.
2. **Move the reference data to a GCS bucket** (Standard storage is roughly 4-5x cheaper
  per GB than a persistent disk) and copy it onto a fresh disk only when you spin up a
  VM. More setup work than what's automated here - worth doing if you settle into a
  regular usage pattern, not built out in these scripts yet.

## Full teardown
When you're done for good (or want to stop the storage cost between long gaps in use):
```bash
gcloud compute instances delete gatk-pipeline-vm --project=<your-project> --zone=us-central1-a
gcloud compute disks delete gatk-pipeline-data --project=<your-project> --zone=us-central1-a
```
This deletes the downloaded references/CADD databases too - you'll redo the one-time
setup next time.
