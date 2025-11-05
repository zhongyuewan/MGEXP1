# Autonomous Reef Monitoring Structures (ARMS) Reveal Human-Induced Biodiversity Shifts in Urban Coastal Ecosystems

#### Zhongyue Wan <sup>1,2,†</sup>, Shelby E. McIlroy <sup>1,2,3,† </sup>, Isis Guibert <sup>1,2</sup>, Wing Yi Haze Chung <sup>1,2,4</sup>, Charlotte Ho <sup>1,2,5</sup>, Alison Corley <sup>1,2</sup>, Emily Chei <sup>1,2</sup>, Inga Conti-Jerpe <sup>6</sup> , Jonathan D. Cybulski <sup>1,2</sup>, Joseph Brennan <sup>1,2</sup>, Ling Fung Matt Chan <sup>1,2</sup>, Philip Thompson <sup>1,2</sup>, Róisín Hayden <sup>1,2</sup>, Shan Yee Joyce Lee <sup>1,2</sup>, Wan Ching Rachel Au <sup>1,2</sup>, Wendy McLeod <sup>1,2</sup>, David M. Baker <sup>1,2,* </sup> 

1 The Swire Institute of Marine Science, The University of Hong Kong, Cape D’Aguilar Road, Shek O, Hong Kong SAR <br>
2 School of Biological Sciences, The University of Hong Kong, Pok Fu Lam, Hong Kong SAR <br>
3 Simon F. S. Li Marine Science Laboratory, School of Life Sciences, The Chinese University of Hong Kong, Shatin, Hong Kong SAR <br>
4 Department of Biology, University of Oxford, Oxford, UK <br>
5 Tree of Life, Wellcome Sanger Institute, Hinxton, Cambridge, UK <br>
6 Lingnan University, Tuen Mun, Hong Kong SAR <br>

`†` These authors contributed equally to the work and should be considered as joint first author. <br>
`*` Corresponding author

## Abstract 

<img align="right" src="2_figure/figure2forshow.png" width=450> 

Urbanization and the global growth of coastal megacities impose an array of stressors on nearby marine ecosystems impacting biodiversity and ecosystem services. Despite extensive research on biodiversity responses to anthropogenic stressors, community level resistance and resilience dynamics – particularly in coastal marine ecosystems – remain poorly understood. In the first application of standardized Autonomous Reef Monitoring Structures (ARMS) as an experimental unit – using a modified common garden experimental design – we examined community-level responses to stress within an urbanized seascape. ARMS that were naturally seeded within marine reserves, accumulating ~1,300 OTUs per unit, across 30 phyla, were transplanted to sites of stress, including domestic sewage discharge and mariculture. We hypothesized that COI metabarcoding of ARMS communities would show that 1) human impacts reduce species richness and alter composition of established communities; and 2) chronic exposure to stressors reduces community resilience following stress removal. We found treatment specific responses and an overall impact of nutrient pollution, particularly inorganic nitrogen, in significantly reducing species richness and altering community structure. Yet, when ARMS were returned to marine reserves, community composition reverted to similarity with ARMS not exposed to stress, accounting for succession. This community-level pattern of low resistance but high resilience highlights the recovery potential of highly connected marine ecosystems and underscores the value of stress mitigation in urbanized seascapes.


For more details, please refer to the [preprint](https://doi.org/10.1101/2025.07.11.664282).   


## Table of Contents

### Supporting Materials 
  1. [Raw sequence](https://doi.org/10.6084/m9.figshare.29481053) 
  2. [Data](3_data)
  3. [Figures](2_figure/250806_merge_figure.pdf)
  4. [Tables](2_figure/250703_merge_table.pdf)
  5. [Supplementary Materials](2_figure/250711_supplementaryMaterials.pdf)

### Sequence processing pipeline 
1. [Import & cutadap](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.1_importAndCutAdapt.sh): import raw sequence data (.fastq) into Qiime artefacts (.qza) and remove PCR adaptors.
2. [Denoise-paired](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.2_denoiseAndPair.sh): remove sequences likely induced by error and merge the reverse/forward reads.
3. [Decontam](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.3_decontam.r): a process to look into the negative control and remove sequences that might have come from sample contamination.
4. [Amino Acid translation](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.4_aaTranslate.r): translate DNA sequence into amino acid and remove sequences with one of the following conditions: 1) any STOP codon, 2) >3 deletion, 3) any frameshift, 4) any insertion.
5. [Cluster all sequences](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.5_clusterReads.sh) by 97% similarity into operational taxonomic units (OTUs) for downstream data analysis.
6. [Taxonomic assignment](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/1.6_taxAssign.sh) with BLAST against two different libraries: 1) McIlroy et al. 2024 & 2) Medori2 (GB260).

### Data Analysis 
1. Environmental data
   - [Heatmap](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.1_eData_heatmap.r) (Figure 1d, Table 1)
   - [MPA east vs west](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.2_eastVSwest.r) (Table S4)
2. Species richness by ARMS 
   - [Merge richness from all three fractions](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.3_combinFractionbyARMS.r) (Table S2)
   - [Environmental data ~ species richness](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.4_eDATAvsRichness.r) (Table 2) 
3. Community composition
   - [PCoA](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.5_PCoA.r) (Figure 2)
   - [Permutational Multivariate Analysis of Variance (adonis2)](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.6_adonis2.r) 
   - [Ternary plot](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.7_ternary.r) (Figure S2)
   - [Horizontal barplot](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.8_sidewayBar.r) (Figure 3a)
   - [Resistance composition](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.9_resistanceComposition.r)(Figure 3b)
4. Succession
   - [Total richness and by phyla](https://github.com/zhongyuewan/MGEXP1/blob/main/1_code/2.10_succession.r) (Figure 4, Table S2)
     
