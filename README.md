# MGEXP1
Pipeline for DNA sequence data &amp; codes for data analysis 


Data processing includes two major parts: 1) converting raw DNA sequence data from metabarcoding into a frequency table and a sequence table & 2) conducting data analysis.

Raw DNA sequence processing follows the following steps:
1. Import raw sequence data (.fastq) into Qiime artefacts (.qza).
2. Remove PCR adaptors.
3. Denoise-paired: remove sequences likely induced by error and merge the reverse/forward reads.
4. Decontam: a process to look into the negative control and remove sequences that might have come from sample contamination.
5. Amino Acid translation: translate DNA sequence into amino acid and remove sequences with any of the following: 1) stop codons, 2) >3 deletion, 3) frameshift, 4) insertion.
6. Cluster all sequences by 97% similarity into operational taxonomic units (OTUs) for downstream data analysis.
7. Taxonomic assignment with BLAST against two different libraries: 1) McIlroy et al. 2024 & Medori2 (GB260).   

Data analysis follows the following steps: 
1. 
