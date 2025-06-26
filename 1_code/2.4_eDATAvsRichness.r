## load lib
library(fmsb)
library(tidyverse)
library(vegan)
library(ggord)
library(pairwiseAdonis)
library(gridExtra)

## set WD 
setwd("/filepath")
rm(list=ls())

# read file i am only getting the middle data, so all sites should have 12 data point (1/month) by the end
metadata.d <- read.csv("../decontamByArms/alphaTable.byArms.edata.csv", header = T)
edata <- read.csv("edata.csv", header = T)
edata.d1 <- edata[,c(31,32,3,4,5,6,26,27,9,15,13,29,22,23,21,25,24,28,19)] %>% filter (phase=="middle")

# check if all sites are evenly assigned a eData
nrow(edata.d1 %>% filter(site=="TPC")) # 30
nrow(edata.d1 %>% filter(site=="CDA")) # 69

# turn all < in to 0
# tunn all NA into na

edata.d2 <- edata.d1
for (i in 7:19){# remove still "<0.1"
edata.d2[,i] <- as.numeric(ifelse(grepl("<", edata.d2[,i]), 0, edata.d2[,i]))
}

# define all rows
site <- c("PC","CC","LM","CDA", "NP","SK","CI","TPC")
phase <- c("initial", "middle","final")

PCrow <- rbind(edata.d2[(edata.d2$Station=="SM10" & edata.d2$Depth =="Bottom Water"),],
               edata.d2[(edata.d2$Station=="SM9" & edata.d2$Depth =="Middle Water"),],
               edata.d2[(edata.d2$Station=="SM11" & edata.d2$Depth =="Middle Water"),])

CCrow <- rbind(edata.d2[(edata.d2$Station=="ST1" & edata.d2$Depth =="Bottom Water"),],
                edata.d2[(edata.d2$Station=="SM7" & edata.d2$Depth =="Middle Water"),],
                edata.d2[(edata.d2$Station=="SM6" & edata.d2$Depth =="Middle Water"),])

LMrow <- rbind(edata.d2[(edata.d2$Station=="SM3" & edata.d2$Depth =="Bottom Water"),],
               edata.d2[(edata.d2$Station=="SM4" & edata.d2$Depth =="Middle Water"),])

CDArow <- rbind(edata.d2[(edata.d2$Station=="SM1" & edata.d2$Depth =="Middle Water"),],
               edata.d2[(edata.d2$Station=="MM8" & edata.d2$Depth =="Surface Water"),])

NProw <- edata.d2[(edata.d2$Station=="MM19" & edata.d2$Depth =="Surface Water"),]
         
SKrow <- edata.d2[(edata.d2$Station=="PM4" & edata.d2$Depth =="Bottom Water"),]

CIrow <- edata.d2[(edata.d2$Station=="TM4" & edata.d2$Depth =="Bottom Water"),]

TPCrow <- edata.d2[(edata.d2$Station=="MM5" & edata.d2$Depth =="Surface Water"),]

# rework PCrow to mean data of the same period from different site into one 
# so all 8 sites have same same data points 

# PCrow
nPCrow <- as.data.frame(matrix(nrow=11, ncol=ncol(PCrow)))
names(nPCrow) <- names(PCrow)
nPCrow[,1:2] <- PCrow[1:11,1:2]
PCrow[,4] <- as.numeric(substr(PCrow[,4],6,7))
nPCrow[,4] <- PCrow[1:11,4]

for (i in 1:11) {
  for (j in 1:13) {
    nPCrow[i,j+6] <- mean((PCrow%>%filter(Dates==nPCrow[i,]$Dates)%>%.[j+6])[,1])
    
  }
  
}


# LMrow
LMrow[,4] <- as.numeric(substr(LMrow[,4],6,7))
nLMrow <- as.data.frame(matrix(nrow=12, ncol=ncol(LMrow)))
names(nLMrow) <- names(LMrow)
nLMrow[,1:2] <- LMrow[1:12,1:2]
nLMrow[,4] <- c(1:12)

for (i in 1:12) {
  for (j in 1:13) {
    nLMrow[i,j+6] <- mean((LMrow%>%filter(Dates==nLMrow[i,]$Dates)%>%.[j+6])[,1])
    
  }
  
}


# CDArow
CDArow[,4] <- as.numeric(substr(CDArow[,4],6,7))
nCDArow <- as.data.frame(matrix(nrow=11, ncol=ncol(CDArow)))
names(nCDArow) <- names(CDArow)
nCDArow[,1:2] <- CDArow[1:11,1:2]
nCDArow[,4] <- c(1,2,4:12)

for (i in 1:11) {
  for (j in 1:13) {
    nCDArow[i,j+6] <- mean((CDArow%>%filter(Dates==nCDArow[i,]$Dates)%>%.[j+6])[,1])
    
  }
  
}

#turn all the row into 
edata.d3 <- rbind(nPCrow,nLMrow,nCDArow,SKrow,CIrow,TPCrow) # 67 rows in total, all sites have 11:12 rows

# write in richness 
metadata <- metadata.d[,1:5]
names(metadata) <- c("ARMS","OTUs","uniOTUs","perUni","site")

edata.d3[,20:21] <- 0
colnames(edata.d3)[20:21] <- c("rich","ARMS")

edata.d4 <- rbind(edata.d3,edata.d3)

for (i in 1:nrow(edata.d3)) {
  for (j in 1:2) {
    
    edata.d4[i+(j-1)*nrow(edata.d3),20] <- metadata%>%filter(site==edata.d3[i,1])%>%.[j,2]
    edata.d4[i+(j-1)*nrow(edata.d3),21] <- metadata%>%filter(site==edata.d3[i,1])%>%.[j,1]
    

  }
  
}

# input treatment and region in edata.d4
edata.d4$treatment <- "MP"
edata.d4$region <- "west"

edata.d4 <- edata.d4 %>% mutate(region= ifelse(site %in% c("TPC","SK","NP","CI"),"east", region))
edata.d4 <- edata.d4 %>% mutate(treatment= ifelse(site %in% c("LM","SK"),"aquaculture", treatment))
edata.d4 <- edata.d4 %>% mutate(treatment= ifelse(site %in% c("CI","PC"),"sewage", treatment))

edata.d4$treatment <- factor(edata.d4$treatment, levels =  c("MP", "aquaculture","sewage"))
edata.d4$site <- factor(edata.d4$site, levels =  c("TPC", "CDA","SK","LM","CI","PC"))

# run a forloop the get all the rich vs eData out 
# rearrange the alpha data to make it easier to loop 
statsDF <- as.data.frame(matrix(0, nrow= 2, ncol= 13))
row.names(statsDF) <- c("pValue","AdjR2")
names(statsDF) <- names(edata.d4)[7:19]
edata.d5 <- edata.d4[1:67,]


# 1,2 it's okay no need to transfer 
# 3: 13 do log(x +1)

for (i in c(1,2)) {
  modelTemp <- lm(rich~ . ,edata.d4[,c("rich",names(statsDF)[i])])
  statsDF[1,i] <- summary(modelTemp)$coefficients[2,4] 
  statsDF[2,i] <- summary(modelTemp)$adj.r.squared
  
}

for (i in c(3:13)) {
  formula_str <- as.formula(paste("log(",names(statsDF)[i], "+1)~ rich"))
  modelTemp <- lm(formula_str, edata.d4)
  statsDF[1,i] <- summary(modelTemp)$coefficients[2,4] 
  statsDF[2,i] <- summary(modelTemp)$adj.r.squared
}  

statsLong <- as.data.frame(t(statsDF))
statsLong$adjPvalueBH <- p.adjust(statsLong$pValue, method = "BH")
statsLong$adjPvalueHOLM <- p.adjust(statsLong$pValue, method = "holm")
statsLong$adjPvalueBonferroni <- p.adjust(statsLong$pValue, method = "bonferroni")


statsLongSorted <- statsLong %>%
  arrange(.[[5]]) 

write.csv(statsLongSorted, "../decontamByArms/noSED/AlphaEdataBreak.csv",row.names = T)
