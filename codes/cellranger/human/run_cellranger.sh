#!/bin/bash

#BSUB -J cellranger_human[1-8]   # Job name and array range
#BSUB -o cellranger_output_%I.out   # Output file
#BSUB -e cellranger_error_%I.err   # Error file
#BSUB -n 20                        # Number of cores
#BSUB -R "rusage[mem=100G]" # Request 100 GB of memory per core, all cores on one machine
#BSUB -q verylong


# Load the Cell Ranger module
module add cellranger/7.1.0

# Set the path to the reference package
reference="/home/i439h/projects/pools/AG_Thongjuea/Dataset/10x_multiomics/references/human/refdata-gex-GRCh38-2024-A"
species="human"

# Set the path to the output directory
output_dir="/home/i439h/projects/hipo_k35/results/cellranger/unforced/$species"

# Define an array of samples
samples=("K35R-BNTVHX" "K35R-BZ8JPB" "K35R-VFHB3C" "K35R-4DZ4CT"  "K35R-7P6VFX" "K35R-sn_ICGC_GBM15"  "K35R-sn_ICGC_GBM71")
sample=${samples[$((LSB_JOBINDEX-1))]}


# Use the LSF job array index to get the right sample
echo $sample
data=/home/i439h/projects/hipo_k35/data/linked_fastqs/$sample
#cd $data

# Run cellranger for the current sample
cellranger count --id="$sample" \
   --transcriptome="$reference" \
  --fastqs="$data" \
  --localcores=20 
#mv $sample $output_dir
