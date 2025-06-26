# What are the files & what do they do 

### [1_edata.csv](1_edata.csv)
Environmental data downloaded directly from [HKEPD](https://cd.epic.epd.gov.hk/EPICRIVER/marine/?lang=en).
The data included all 13 environmental parameters used in the analysis. 

### [2_sample-metadata.tsv](2_sample-metadata.tsv)
Raw metadata including the sequence path for bioinformatics.

### [3_sample-metadata.byArms.tsv](3_sample-metadata.byArms.tsv)
Metadata after merging the fractions. 

### [4_feaTabClean.csv](4_feaTabClean.csv)
Frequency table with all the sequences before clustering by 97% similarity by fractions. 

### [5_freqTableCleanByArms.csv](5_freqTableCleanByArms.csv)
Frequency table after clean up and clustered by ARMS.

### [6_dna-sequences.fasta](6_dna-sequences.fasta)
All raw sequences before decontamination and clustering.

### [7_TaxAsn_shelbyOmidori.csv](7_TaxAsn_shelbyOmidori.csv)
Taxonomically assigned sequence. Results came from alignment against two databases (one locally curated and Midori2).
