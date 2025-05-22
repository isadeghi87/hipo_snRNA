#!/bin/bash

#BSUB -J cellranger_array[1-2]   # Job name and array range
#BSUB -o cellranger_output_%I.out   # Output file
#BSUB -e cellranger_error_%I.err   # Error file
#BSUB -n 20                        # Number of cores
#BSUB -R "rusage[mem=150G]" # Request 150 GB of memory per core, all cores on one machine
#BSUB -q verylong


# Load the Cell Ranger module
module add cellranger/7.1.0

# Set the path to the reference package
reference="/home/i439h/projects/pools/AG_Thongjuea/Dataset/10x_multiomics/references/human_mouse/refdata-gex-GRCh38_and_GRCm39-2024-A"
species='PDX'

# Set the path to the output directory
output_dir="/home/i439h/projects/hipo_k35/results/cellranger/unforced/$species"

# Define an array of samples
#samples=("K35R-BR284Y"  "K35R-Tr041"  "K35R-Tr068"  "K35R-Tr082"  "K35R-Tr083"  "K35R-XN9L9Y"  "K35R-sn_Astro110fh"  "K35R-sn_SU_DIPG_13")
samples=("K35R-sn_Astro110fh")
sample=${samples[$((LSB_JOBINDEX-1))]}


# Use the LSF job array index to get the right sample
echo $sample
data="/home/i439h/projects/hipo_k35/data/linked_fastqs/$sample"
#cd $data

# Run cellranger for the current sample
cellranger count --id="$sample" \
   --transcriptome="$reference" \
  --fastqs="$data" \
  --localcores=20 
#mv $sample $output_dir
