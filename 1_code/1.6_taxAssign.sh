#!/bin/bash

#SBATCH --job-name=taxassigBLAST                  # Job name
#SBATCH --partition=amd                           # Specific Partition (intel/amd)
#SBATCH --nodes=1                                 # number of nods
#SBATCH --qos=normal                              # Specific QOS (debug, normal)
#SBATCH --time=12:00:00                           # Wall time limit (days-hrs:min:sec)
#SBATCH --cpus-per-task=2                         # Number of cpu/task 
#SBATCH --mem=32G                                 # Request total amount of RAM
#SBATCH --output=%x.out.%j                        # Standard output file
#SBATCH --error=%x.err.%j                         # Standard error file


###### THE FOLLOWING IS MY SHIT
# load environment and qiime2 
module load blast-plus/2.13.0

# converting db 
echo "coverting DB `date "+%Y/%m/%d -- %H:%M:%S"`"
makeblastdb -in 10-taxAssign/MIDORI2_LONGEST_NUC_SP_GB260_CO1_SINTAX.fasta -dbtype nucl -out 10-taxAssign/MIDORI_DB_BLAST/taxonomy_db
makeblastdb -in 10-taxAssign/shelby_SINTAX.fasta -dbtype nucl -out 10-taxAssign/Shelby_DB_BLAST/taxonomy_db


# tax assignment 
echo "Tax assignment `date "+%Y/%m/%d -- %H:%M:%S"`"
blastn -query 10-taxAssign/dna-sequences.fasta -db 10-taxAssign/MIDORI_DB_BLAST/taxonomy_db -out 10-taxAssign/blast_midori.txt -outfmt 6
blastn -query 10-taxAssign/dna-sequences.fasta -db 10-taxAssign/Shelby_DB_BLAST/taxonomy_db -out 10-taxAssign/blast_shelby.txt -outfmt 6


echo "finished by `date "+%Y/%m/%d -- %H:%M:%S"`"

