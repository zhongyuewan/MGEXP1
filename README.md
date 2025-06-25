# Autonomous Reef Monitoring Structures (ARMS) Reveal Human-Induced Biodiversity Shifts in Urban Coastal Ecosystems

Author list 
Abstract here

# Table of Contents

- raw sequence data from [xxx](link)
- metadata


Data processing includes two major parts: 1) converting raw DNA sequence data from metabarcoding into a frequency table and a sequence table & 2) conducting data analysis.

Raw DNA sequence processing follows the following steps:
1. [Import raw sequence data (.fastq) into Qiime artefacts (.qza) and remove PCR adaptors.](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.1_importAndCutAdapt.sh)
2. [Denoise-paired: remove sequences likely induced by error and merge the reverse/forward reads.](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.2_denoiseAndPair.sh)
3. [Decontam: a process to look into the negative control and remove sequences that might have come from sample contamination.](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.3_decontam.r)
4. [Amino Acid translation: translate DNA sequence into amino acid and remove sequences with any of the following: 1) stop codons, 2) >3 deletion, 3) frameshift, 4) insertion.]
5. [Cluster all sequences by 97% similarity into operational taxonomic units (OTUs) for downstream data analysis.]
6. [Taxonomic assignment with BLAST against two different libraries: 1) McIlroy et al. 2024 & Medori2 (GB260).]

Data analysis follows the following steps: 
1. 
