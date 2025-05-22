# HIPO_snRNA

**Single-Nucleus RNA-Seq Analysis repo for pediatric high-grade glioma cell composition from allograft**

## Overview

This repository contains a complete end-to-end workflow for processing and analyzing single-nucleus RNA sequencing (snRNA-seq) data from allograft. The pipeline is designed to:

- Perform raw data QC and filtering
- Align nuclei-level reads and generate expression matrices
- Conduct downstream analysis including normalization, clustering, and cell-type annotation
- Perform differential expression, sub-clustering, and trajectory inference
- Visualize results via interactive notebooks and publication-quality plots

## Key Features

- **Automated preprocessing**: demultiplexing and alignment with Cell Ranger
- **Quality control**: per-nucleus metrics, doublet detection, and filtering
- **Comprehensive analysis**: normalization, dimensionality reduction, clustering, and annotation using Seurat
- **Advanced downstream**: marker gene identification, pseudotime/trajectory analysis with Monocle3 or Slingshot
- **Reproducible reports**: RMarkdown and Jupyter notebooks for each analysis step

## Requirements

- **Cell Ranger** (10x Genomics) v6.0 or higher
- **Conda** (miniconda or Anaconda)
- **Python** ≥ 3.8
  - `scanpy`, `anndata`, `numpy`, `pandas`, `scikit-learn`, `matplotlib`
- **R** ≥ 4.1
  - `Seurat` (v4), `tidyverse`, `DoubletFinder`, `Monocle3` or `Slingshot`, `patchwork`
- **RMarkdown** for generating HTML/PDF reports

## Installation

1. **Clone the repository**

   ```bash
   git clone https://github.com/isadeghi87/hipo_snRNA.git
   cd hipo_snRNA
   ```

2. **Set up Conda environments**

   ```bash
   # Create environments
   conda create -n hipo_snrna_py python=3.9 scanpy anndata numpy pandas scikit-learn matplotlib
   conda create -n hipo_snrna_r r-base=4.1 seaborn
   
   # Activate and install R packages
   conda activate hipo_snrna_r
   Rscript -e "install.packages(c('Seurat', 'tidyverse', 'patchwork', 'DoubletFinder', 'Monocle3', 'Slingshot'))"
   ```

3. **Install Cell Ranger**

   Follow the official instructions at https://support.10xgenomics.com.

## Directory Structure

```text
├── raw_data/                # FASTQ files organized by sample
├── preprocessing/           # Scripts for demux, alignment, and count matrix generation
│   └── cellranger/          # Cell Ranger config and wrapper scripts
├── qc/                      # Quality control scripts and per-nucleus metrics
├── analysis/                # R and Python scripts for downstream processing
│   ├── seurat_pipeline.R    # Normalization, clustering, annotation
│   ├── doublet_detection.R  # Identify and remove doublets
│   └── trajectory_analysis/ # Monocle3/Slingshot workflows
├── notebooks/               # Jupyter notebooks for interactive exploration
├── reports/                 # RMarkdown templates and generated HTML/PDF
└── README.md                # This file
```

## Usage

### 1. Preprocessing

Demultiplex and align raw FASTQs:

```bash
bash preprocessing/cellranger/run_cellranger.sh \
  --fastqs raw_data/Sample1/ \
  --id Sample1 \
  --transcriptome /path/to/refdata-cellranger
```

### 2. Quality Control

Generate QC metrics and filter nuclei:

```bash
Rscript qc/compute_qc_metrics.R --input analysis/Sample1/raw_matrix.rds --output qc/Sample1_qc_metrics.csv
Rscript qc/filter_nuclei.R --metrics qc/Sample1_qc_metrics.csv --output analysis/Sample1/filtered_matrix.rds
```

### 3. Downstream Analysis

Run the Seurat pipeline for clustering and annotation:

```bash
Rscript analysis/seurat_pipeline.R \
  --input analysis/Sample1/filtered_matrix.rds \
  --output analysis/Sample1/seurat_obj.rds
```

Perform doublet detection:

```bash
Rscript analysis/doublet_detection.R --seurat_obj analysis/Sample1/seurat_obj.rds
```

### 4. Trajectory Inference

Execute trajectory workflows:

```bash
Rscript analysis/trajectory_analysis/run_monocle3.R \
  --seurat_obj analysis/Sample1/seurat_obj.rds \
  --output analysis/Sample1/monocle3_trajectory.rds
```

### 5. Interactive Exploration

Launch Jupyter notebooks:

```bash
conda activate hipo_snrna_py
jupyter notebook notebooks/
```

## Outputs

- **Count matrices**: raw and filtered `.rds`/`.h5ad` files
- **QC reports**: CSV summaries and HTML visualizations
- **Seurat objects**: annotated cell clusters and metadata
- **Trajectory objects**: pseudotime assignments and lineage trees
- **Publication figures**: UMAP/t-SNE plots, violin and feature plots, trajectory graphs

## Contributing

Contributions welcome! To contribute:

1. Fork this repository
2. Create a topic branch: `git checkout -b feature/your-feature`
3. Commit your changes: `git commit -m "Add feature"
4. Push and open a Pull Request

Please follow existing code style and include tests or example notebooks where appropriate.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
