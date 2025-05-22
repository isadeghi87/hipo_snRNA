#!/bin/bash

#BSUB -J cellranger_mouse[1-7]   # Job name and array range
#BSUB -o cellranger_output_%I.out   # Output file
#BSUB -e cellranger_error_%I.err   # Error file
#BSUB -n 20                        # Number of cores
#BSUB -R "rusage[mem=150G]" # Request 100 GB of memory per core, all cores on one machine
#BSUB -q verylong


# Load the Cell Ranger module
module add cellranger/7.1.0

# Set the path to the reference package
reference="/home/i439h/projects/pools/AG_Thongjuea/Dataset/10x_multiomics/references/mouse/GRCm39/refdata-gex-GRCm39-2024-A/"
species='mouse'

# Set the path to the output directory
output_dir="/home/i439h/projects/hipo_k35/results/cellranger/unforced/$species"

# Define an array of samples
#samples=("K35R-3Q7AGR"  "K35R-A39LMU"  "K35R-EP019"  "K35R-RERPK9"  "K35R-Tr085"  "K35R-Tr086"  "K35R-Tr089"  "K35R-Tr091"  "K35R-Tr103"  "K35R-Tr107"  "K35R-Tr110")
samples=("K35R-CGHTYL" "K35R-EGNNNG" "K35R-RTJDH5" "K35R-D2BCNW" "K35R-MLGFBW" "K35R-TRL883" "K35R-8PAV6H")
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
