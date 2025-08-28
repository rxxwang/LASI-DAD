mkdir -p log
snakemake \
   --keep-going \
   --jobs 500 \
   --max-jobs-per-second 5 \
   --latency-wait 60 \
   --cluster-config cluster.yaml  \
   --cluster "sbatch \
              --output={cluster.log}_%j.out \
              --error={cluster.log}_%j.err \
              --job-name={cluster.name} \
              --time={cluster.time}  \
              --cpus-per-task={cluster.cpus}  \
              --mem={cluster.mem}"
