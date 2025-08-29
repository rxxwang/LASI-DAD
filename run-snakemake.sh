snakemake \
   --keep-going \
   --jobs 500 \
   --max-jobs-per-second 5 \
   --latency-wait 60 \
   --executor slurm \
   --default-resources \
   --slurm-logdir /net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/snake-log/ \
   --slurm-keep-successful-logs
