#!/bin/bash

#SBATCH --job-name=cluster97                    # Job name
#SBATCH --partition=amd                           # Specific Partition (intel/amd)
#SBATCH --nodes=1                                 # number of nods
#SBATCH --qos=normal                              # Specific QOS (debug, normal)
#SBATCH --time=12:00:00                           # Wall time limit (days-hrs:min:sec)
#SBATCH --cpus-per-task=1                         # Number of cpu/task 
#SBATCH --mem=32G                                 # Request total amount of RAM
#SBATCH --output=%x.out.%j                        # Standard output file
#SBATCH --error=%x.err.%j                         # Standard error file


###### THE FOLLOWING IS MY SHIT
# load environment and qiime2 
module load QIIME2/2021.11

# do this after AA translation 
#convert freTable and rep=seqs back to qza 
echo "0/3 convert things to .biom"
cd /lustre1/g/sbs_dmb/wilson/20240529_MG_exp1_all/8.6-AAtranslatED
biom convert -i feaTabClean.txt -o feaTabClean.biom --table-type="OTU table" --to-json
cd ..

echo "1/3 truning the decontaminate data back to .qza"
qiime tools import \
  --input-path 8.6-AAtranslatED/feaTabClean.biom \
  --output-path 8.6-AAtranslatED/feaTabClean.qza \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format 

qiime tools import \
  --input-path 8.6-AAtranslatED/repSeqClean.fasta \
  --output-path 8.6-AAtranslatED/repSeqClean.qza \
  --type 'FeatureData[Sequence]'

#clustering of ASVs into OTUs at 97% - this differs among genes
echo "2/3 cluster to 97%"
qiime vsearch cluster-features-de-novo \
  --i-table 8.6-AAtranslatED/feaTabClean.qza \
  --i-sequences 8.6-AAtranslatED/repSeqClean.qza \
  --p-perc-identity 0.97 \
  --o-clustered-table table-dn-97.qza \
  --o-clustered-sequences rep-seqs-dn-97.qza
  
