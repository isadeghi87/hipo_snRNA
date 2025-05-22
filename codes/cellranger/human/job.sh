bsub -R "rusage[mem=150G]" -n 20 -q verylong < run_cellranger.sh
