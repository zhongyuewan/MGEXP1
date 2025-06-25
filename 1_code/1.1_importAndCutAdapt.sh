#!/bin/bash

#SBATCH --job-name=importAndCutAdapt              # Job name
#SBATCH --partition=amd                           # Specific Partition (intel/amd)
#SBATCH --qos=normal                              # Specific QOS (debug, normal)
#SBATCH --time=12:00:00                           # Wall time limit (days-hrs:min:sec)
#SBATCH --ntasks=4                                # Number of tasks
#SBATCH --mem=32G                                 # Request total amount of RAM
#SBATCH --output=%x.out.%j                        # Standard output file
#SBATCH --error=%x.err.%j                         # Standard error file


###### THE FOLLOWING IS MY SHIT
# load environment and qiime2 
source 1-scripts/0-module_load

# 1. import all the sequence data from fastq
echo "import file to qiime"

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path 0-metadata/sample-metadata.tsv \
  --output-path 4-imported/paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

echo "summary file"
qiime demux summarize \
  --i-data 4-imported/paired-end-demux.qza \
  --o-visualization 4-imported/paired-end-demux-summ.qzv

# 2. cutadapter
# note: 
# TAIACYTCIGGRTGICCRAARAAYCA is replaced by 
# TANACYTCNGGRTGNCCRAARAAYCA because I is inosine which is not in IUPC
# Inosine binds with A/C/U(T), so i'm changing it to N. 

echo "cutadapt "

qiime cutadapt trim-paired \
  --p-cores 20 \
  --i-demultiplexed-sequences 4-imported/paired-end-demux.qza \
  --p-front-f GGWACWGGWTGAACWGTWTAYCCYCC\
  --p-front-r TANACYTCNGGRTGNCCRAARAAYCA\
  --p-error-rate 0 \
  --o-trimmed-sequences 5-cutadapt/trimmed-seqs.qza \
  --verbose

echo "visulize data"
qiime demux summarize \
--i-data 5-cutadapt/trimmed-seqs.qza \
--o-visualization 5-cutadapt/trimmed-seqs.qzv
