# Autonomous Reef Monitoring Structures (ARMS) Reveal Human-Induced Biodiversity Shifts in Urban Coastal Ecosystems

#### Zhongyue Wan 1,2 , Isis Guibert 1,2 , Wing Yi Haze Chung 1,2,4 , Charlotte Ho 1,2,5 , Alison Corley 1,2 , Emily Chei 1,2 , Inga Conti-Jerpe 1,2,6 , Jonathan D. Cybulski 1,2 , Joseph Brennan 1,2 , Ling Fung Matt Chan 1,2 , Philip Thompson 1,2 , Róisín Hayden 1,2 , Shan Yee Joyce Lee 1,2 , Wan Ching Rachel Au 1,2 , Wendy McLeod 1,2 , David M. Baker 1,2,* &amp; Shelby E. McIlroy 1,2,3,*

1 The Swire Institute of Marine Science, The University of Hong Kong, Cape D’Aguilar Road, Shek O, Hong Kong SAR <br>
2 School of Biological Sciences, The University of Hong Kong, Pok Fu Lam, Hong Kong SAR <br>
3 Simon F. S. Li Marine Science Laboratory, School of Life Sciences, The Chinese University of Hong Kong, Shatin, Hong Kong SAR <br>
4 University of Oxford, Oxford, UK <br>
5 Tree of Life, Wellcome Sanger Institute, Hinxton, Cambridge, UK <br>
6 Lingnan University, Tuen Mun, Hong Kong SAR <br>

`*` Corresponding author

Biodiversity thrives in coastal marine habitats which host foundational species such as corals, mangroves, and seagrasses. However, coastal development and the growth of megacities along shorelines impose an array of stressors on the marine environment. These stressors inevitably impact biodiversity which dictates ecosystem functions and services. Despite extensive research on biodiversity responses to anthropogenic stressors, phylum-specific resistance and resilience dynamics – particularly in coastal marine ecosystems – remain poorly understood. Considering the global scale of coastal development, it is imperative to develop a more comprehensive understanding of how biodiversity, in terms of richness and community composition, is influenced by various anthropogenic stressors. Here, we present the first application of standardized Autonomous Reef Monitoring Structures (ARMS) as an experimental unit - using a common garden experimental design - to examine community responses to stress within an urbanized seascape. ARMS were seeded within a marine reserves for one year and then transplanted to sites of stress, including domestic sewage, and mariculture. We hypothesized that 1) human impacts alter the richness and composition of established communities; and 2) the intensity of these impacts influences community resistance and resilience to stress. We found that although community resistance to stress was low - particularly under high-nutrient stress - resilience was high, suggesting that urbanized seascapes have high recovery potential when stress is mitigated.

For more details, please refer to the [preprint_pending](link)   


## Table of Contents

### Raw sequence data from [xxx](link)

### Supporting Materials 
  1. [Metadata](link)
  2. [Frequency table](link)
  3. [Sequence table](link)
  4. [Taxonomic assignment](link)
  5. [Figures](link)
  6. [Tables](link)
  7. [Supplementary Materials](link)

### All codes for sequence processing and data analysis. Specifically:
- Raw DNA sequence processing follows the following steps:
1. [Import raw sequence data (.fastq) into Qiime artefacts (.qza) and remove PCR adaptors.](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.1_importAndCutAdapt.sh)
2. [Denoise-paired](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.2_denoiseAndPair.sh): remove sequences likely induced by error and merge the reverse/forward reads.
3. [Decontam](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.3_decontam.r): a process to look into the negative control and remove sequences that might have come from sample contamination.
4. [Amino Acid translation](link): translate DNA sequence into amino acid and remove sequences with any of the following: 1) stop codons, 2) >3 deletion, 3) frameshift, 4) insertion.
5. [Cluster all sequences](link) by 97% similarity into operational taxonomic units (OTUs) for downstream data analysis.
6. [Taxonomic assignment](link) with BLAST against two different libraries: 1) McIlroy et al. 2024 & Medori2 (GB260).

- Data analysis follows the following steps: 
1. 
