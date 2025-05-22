bsub -R "rusage[mem=100G]" -n 10 -q verylong < run_cellranger.sh
