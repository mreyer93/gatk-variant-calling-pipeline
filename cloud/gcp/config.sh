#!/usr/bin/env bash
# Shared settings for the cloud/gcp/*.sh scripts. Edit these to taste, then source this
# file from the other scripts (they do this automatically).

# Your GCP project ID (`gcloud config get-value project` if you're not sure, or
# `gcloud projects list` to see what you have access to).
PROJECT_ID="${PROJECT_ID:-your-gcp-project-id}"

# Zone to create resources in. us-central1 tends to have good spot availability/pricing.
ZONE="${ZONE:-us-central1-a}"

VM_NAME="${VM_NAME:-gatk-pipeline-vm}"
DISK_NAME="${DISK_NAME:-gatk-pipeline-data}"

# 16 vCPU / 64G RAM - a reasonable default for parallel Mutect2/HaplotypeCaller jobs.
# Adjust based on how many samples/jobs you run concurrently (see --jobs/--cores when
# launching snakemake). Bigger machine = faster wall-clock, same total cost either way
# for a fixed amount of work, since spot billing is per-second.
MACHINE_TYPE="${MACHINE_TYPE:-n2-standard-16}"

# Data disk size in GB. Needs to fit: CADD databases (400G+), the GRCh38/dbSNP/gnomAD/
# Funcotator bundle (~70G), plus working space for pipeline outputs. 1000G gives some
# headroom; resize upward later with `gcloud compute disks resize` if you need more.
DISK_SIZE_GB="${DISK_SIZE_GB:-1000}"
DISK_TYPE="${DISK_TYPE:-pd-balanced}"

BOOT_DISK_SIZE_GB="${BOOT_DISK_SIZE_GB:-50}"
IMAGE_FAMILY="${IMAGE_FAMILY:-debian-12}"
IMAGE_PROJECT="${IMAGE_PROJECT:-debian-cloud}"
