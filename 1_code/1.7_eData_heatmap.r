## load lib
library(tidyverse)
library(factoextra)
library(RColorBrewer)
library(gridExtra)
library(pheatmap)
library(MetBrewer)

## set WD 
# mac
setwd("filepath")
rm(list=ls())

# read file 
edata <- read.csv("edata.csv", header = T)

names(edata) # check what i need 
# General -> Dissolved.Oxygen..mg[26], Dissolved.Oxygen...saturation[27]
# General -> Salinity[18], pH[20]
# sedimentation ->  Turbidity [9], Suspended.Solids[15]
# aquaculture -> Total.Inorganic.Nitrogen[13],Ammonia.Nitrogen[29],
# aquaculture -> Nitrite.Nitrogen[22],Nitrate.Nitrogen[23],Orthophosphate.Phosphorus[21],
# sewage -> E..coli..cfu[25], Faecal.Coliforms[24]
# sewage -> Chlorophyll.a..[28], Phaeo.pigments [19]

# also check Silica for bacillariophyta 

edata.d1 <- edata[,c(31,32,3,4,5,6,13,29,22,23,21,25,24,28,19,26,27,15,9)]

siteP.1 <- read.csv("sitePofile.csv")
siteP <- siteP.1%>% filter(site != "NP" & site != "CC")

# filter dates 
edata.d2 <- edata.d1 %>%
  filter (phase == "middle")
  
# filter depth, 
edata.d3 <- edata.d2 %>%
  mutate(mergecell=paste0(Station,Depth)) %>% 
  filter(mergecell %in% siteP$mergecell)

# turn anything "<" into 0. because the <xx.xx just means the number is too 
# small to be picked by the machine so mathmetically i will treat it as 0
edata.d4 <- edata.d3
for (i in 7:19){# remove "<"
  edata.d4[,i] <- as.numeric(ifelse(grepl("<", edata.d3[,i]), 0, edata.d3[,i]))
}

# mutate a unique id for each sampling point 
edata.d5 <- mutate(edata.d4, ID = paste0(site, row_number()))
row.names(edata.d5) <- edata.d5$ID

# build a new df to store the mean and sd 
edata.d6 <-   data.frame(matrix(0, nrow=8, ncol=13))
row.names(edata.d6) <-  c(unique(edata.d5$site), "mean", "sd")
names(edata.d6) <- names(edata.d5)[7:19]

# make a function to mean all the things 
meanEdata <- function (df){
  for (i in 1:6){
for (j in 1:13){
  edata.d6[i,j] <- mean(df[df$site == row.names(edata.d6[i,]),][,j+6])
}
  }
  return(edata.d6)
}


edata.d6 <- meanEdata(edata.d5)

# calculate group mean and sd
for (i in 1:13){ #mean
edata.d6[7,i] <- mean(edata.d6[1:6,i])

}

for (i in 1:13){ #sd
  edata.d6[8,i] <- sd(edata.d6[1:6,i])
  
}


# turn this table into -mena / sd 
edata.d7 <- edata.d6
for (i in 1:6){ #how many of sd
  for (j in 1:13){
  edata.d7[i,j] <- (edata.d6[i,j] - edata.d6[7,j])/edata.d6[8,j]
}
}



# add treatment/water  
edata.d7$water <- "west"
edata.d7$water[row.names(edata.d7) == "NP" | row.names(edata.d7) == "TPC" |
                 row.names(edata.d7) == "CI" |row.names(edata.d7) == "SK"] <- "east"

edata.d7$treatment <- "baseline"
edata.d7$treatment[row.names(edata.d7) == "PC"|row.names(edata.d7) == "CI"] <- "eutrophication"
edata.d7$treatment[row.names(edata.d7) == "LM"|row.names(edata.d7) == "SK"] <- "fishfarm"

edata.d7$water[7:8] <- "meansd"
edata.d7$treatment[7:8] <- "meansd"

edata.d7 <- mutate(edata.d7, treatmentWater = paste0(water,treatment))

# make clear df to plot figure
edata.d8 <- edata.d7[1:6,] # remove mean/sd
edata.fi <- edata.d8
edata.fi[,10:11] <- -edata.fi[,10:11] # reverse 1,2 which is o2, salinity, pH 

# reorder them by east/west & impact 
customOrder <- c("TPC","SK","CI","CDA","LM","PC")
edata.plot <- edata.fi[match(customOrder,rownames(edata.fi)),]
edata.plot$Mix_impact <- rowMeans(edata.plot[,1:13])



# heat map 
# show color 
myColor <- colorRampPalette(c("#0072B2", "#F0E442", "#D55E00"))(80)


# plot
# add new names 
names(edata.plot)[1:13] <- c("Total Inorganic Nitrogen (mg/l)","Ammonia Nitrogen (mg/l)", "Nitrite Nitrogen (mg/l)", "Nitrate Nitrogen (mg/l)", 
                             "Orthophosphate Phosphorus (mg/l)", "E. coli (cfu/100ml)", "Faecal Coliforms (cfu/100ml)", "Chlorophyll a (µg/L)", "Phaeo pigments (µg/L)",
                             "Dissolved Oxygen  (mg/L)", "Dissolved Oxygen saturation (% of saturation)", "Suspended Solids (mg/L)", "Turbidity (NTU)")

plot1 <- pheatmap(as.matrix(t(edata.plot[,c(1:13)])), cluster_cols = F, cluster_rows = F, 
                 scale="row", col = myColor,fontsize = 20) 


plot1.2<- heatmap(as.matrix(t(edata.plot[,c(7:10,12,13,13)])), Colv = T, Rowv = NA, 
                  scale="row", col = myColor) 



# heatmap side by site
plot.fin.west <- pheatmap(as.matrix(t(edata.plot[4:6,c(1:13)])), cluster_rows = F, cluster_cols = F, 
                        scale="row", col = myColor, fontsize_row = 20) 


plot.fin.east <- pheatmap(as.matrix(t(edata.plot[1:3,c(1:13)])), cluster_rows = F, cluster_cols = F, 
                         scale="row", col = myColor, fontsize_row = 20) 




##### write stuff out 
write.csv(edata.plot,"edata.plot.csv", row.names = T)
write.csv(edata.d6,"edata.plotBeforeNormal.csv", row.names = T)

