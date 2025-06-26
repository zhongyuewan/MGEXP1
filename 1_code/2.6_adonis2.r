## set WD 

#mac 
setwd("filepath")

rm(list=ls())

## load lib
library(dplyr)
library(ggplot2)
library(Biostrings)
library(factoextra)
library(vegan)
library(pairwiseAdonis)


## get data
# it could be jaccard or BC.H, just load the differernt diversity matrix 
# BC.H here is not really BC h, it is BC on square root transfered count data
freqTable <- read.csv('decontamByArms/noSED/freqTable.csv', header = TRUE, row.names = 1)
metadata <- read.csv('decontamByArms/noSED/metadata.csv', row.names = 1)

## fix table names and stuffs 
names(freqTable) <- gsub("X","", names(freqTable))



# check out sharing OTUs by time 
iARMS <- row.names(metadata[(metadata$phase=="initial"),])
mARMS <- row.names(metadata[(metadata$phase=="middle"),])
fARMS <- row.names(metadata[(metadata$phase=="final"),])
baseARMS <- row.names(metadata[(metadata$phase=="initial") | (metadata$phase=="middle") | (metadata$phase=="final"),])
base12 <- row.names(metadata[(metadata$phase=="initial") | (metadata$phase=="middle"),])
base23 <- row.names(metadata[(metadata$phase=="middle") | (metadata$phase=="final"),])
base13 <- row.names(metadata[(metadata$phase=="initial") | (metadata$phase=="final"),])


RTARMS <- row.names(metadata[(metadata$phase=="middle" | metadata$phase=="resistance"),])
RLARMS <- row.names(metadata[(metadata$phase=="final" | metadata$phase=="resilience"),])

wARMS <- row.names(metadata[(metadata$waterAct=="west"),])
eARMS <- row.names(metadata[(metadata$waterAct=="east"),])

wRLARMS <- RLARMS[c(RLARMS %in% wARMS)]
eRLARMS <- RLARMS[c(RLARMS %in% eARMS)]

wRTARMS <- RTARMS[c(RTARMS %in% wARMS)]
eRTARMS <- RTARMS[c(RTARMS %in% eARMS)]

wiARMS <- iARMS[c(iARMS %in% wARMS)]
eiARMS <- iARMS[c(iARMS %in% eARMS)]

wmARMS <- mARMS[c(mARMS %in% wARMS)]
emARMS <- mARMS[c(mARMS %in% eARMS)]

wfARMS <- fARMS[c(fARMS %in% wARMS)]
efARMS <- fARMS[c(fARMS %in% eARMS)]

# Calculate distance matrix 
distMBCsr <-as.matrix(vegdist(t(sqrt(freqTable)), method = "bray")) 

# adonis
# succession, compare the baseline arms over 3 years
# treatment p<0.05 succession is real
# water on treatment p=0.002, east/west impacts on how succession pattern changes 
set.seed(123)
model1 <- adonis2(distMBCsr[baseARMS,]~treatment*waterAct, 
                    data=metadata[baseARMS,],	by = "terms")
model1

set.seed(123)
model1.ph <- pairwise.adonis2(distMBCsr[baseARMS,]~treatment*waterAct, 
                                     data=metadata[baseARMS,],	by = "terms")
model1.ph


# what about resistance treatment
# perfect results 
# water is also a significant impact but no interation from water on treatment 
set.seed(123)
model2 <- adonis2(t(distMBCsr)[RTARMS,]~treatment*waterAct, 
                  data=metadata[RTARMS,],	by = "terms") # water interacts with treatment shifts BCD
set.seed(123)
model2.ph <- pairwise.adonis2(t(distMBCsr)[RTARMS,] ~ treatment*waterAct, 
                              data = metadata[RTARMS,], 
                              p.adjust.method = "BH", permutations=1000,	by = "terms")




# what about resilience treatment? 
# treatment has no imapcts on community structure meaning they all recover back -> high resiliance 
# water however has significant impact on community structure
set.seed(123)
model3 <- adonis2(t(distMBCsr)[RLARMS,]~treatment+waterAct, 
                  data=metadata[RLARMS,],	by = "terms") # water interacts with treatment shifts BCD
set.seed(123)
model3.ph <- pairwise.adonis2(t(distMBCsr)[RLARMS,] ~ treatment*waterAct, 
                              data = metadata[RLARMS,], 
                              p.adjust.method = "BH",	by = "terms")


# what about mbaseline vs all final 
# treatment has no imapcts on community structure meaning they all recover back -> high resiliance 
# water however has significant impact on community structure
metadata.m3 <- metadata[c(mARMS,RLARMS),]
metadata.m3[metadata.m3$phase=="final",]$phase <- "resilience"


set.seed(123)
model3.mid <- adonis2(t(distMBCsr)[c(mARMS,RLARMS),]~phase*waterAct, 
                  data=metadata.m3,	by = "terms") # water interacts with treatment shifts BCD



