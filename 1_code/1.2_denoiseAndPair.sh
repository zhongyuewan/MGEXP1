#!/bin/bash

#SBATCH --job-name=denoiseAndPair     # Job name
#SBATCH --partition=amd               # Specific Partition (intel/amd)
#SBATCH --qos=normal                  # Specific QOS (debug, normal)
#SBATCH --time=12:00:00               # Wall time limit (days-hrs:min:sec)
#SBATCH --ntasks=16                   # Number of tasks
#SBATCH --mem=64G                     # Request total amount of RAM
#SBATCH --output=%x.out.%j            # Standard output file
#SBATCH --error=%x.err.%j             # Standard error file


###### THE FOLLOWING IS MY SHIT
# load environment and qiime2 
# i will use 4 cpu to speed it up
# it took 1 hour for it
# so in the future, i might do 8 cpu and it will take 30min

## new update using 16 cpu and 64G ram 
## see how long does it take

source 1-scripts/0-module_load
echo "denoise and merge"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs 5-cutadapt/trimmed-seqs.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 220 \
  --p-n-threads 16 \
  --o-table 6-denoiseAndMerged/table.qza \
  --o-representative-sequences 	6-denoiseAndMerged/rep-seqs.qza \
  --o-denoising-stats 	6-denoiseAndMerged/denoising-stats.qza

qiime metadata tabulate \
  --m-input-file 	6-denoiseAndMerged/denoising-stats.qza \
  --o-visualization 	6-denoiseAndMerged/denoising-stats.qzv

qiime feature-table summarize \
  --i-table 6-denoiseAndMerged/table.qza \
  --o-visualization 6-denoiseAndMerged/table.qzv \
  --m-sample-metadata-file 0-metadata/sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data 6-denoiseAndMerged/rep-seqs.qza \
  --o-visualization 6-denoiseAndMerged/rep-seqs.qzv
