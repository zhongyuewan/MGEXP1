#########################
# PREPARE R ENVIRONMENT #
#########################
library(Biostrings)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(forcats)
library(tidyr)
library(cowplot)



## set WD 
setwd("/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#1 biodiveristy/#2 CRF EXP1/#2 dataAnalysis_all/1_data")
rm(list=ls())

## get data
freqTable <- read.csv('decontamByArms/noSED/freqTable.csv', header = TRUE, row.names = 1)
metadata <- read.csv('decontamByArms/noSED/metadata.csv', row.names = 1)
sequenceTable <- readDNAStringSet('decontamByArms/noSED/dna-sequences97.fasta')
SOM.d1 <- read.csv('taxAssign/blast/TaxAsn_shelbyOmidori.csv') 

SOM <- SOM.d1 %>% filter(id %in% row.names(freqTable))
nrow(SOM.d1) - nrow(SOM) # 408 sedimentation rows 

#add a $year to show 1year, 2years, 3years arms 
metadata$years <- "3 years"
metadata[metadata$phase=="initial",]$years <- "1 year"
metadata[metadata$phase=="middle" | metadata$phase=="resistance" ,]$years <- "2 years"

#touch up table 
names(freqTable) <- gsub("X","", names(freqTable))


# check out sharing OTUs by time 
iARMS <- row.names(metadata[(metadata$phase=="initial"),])
mARMS <- row.names(metadata[(metadata$phase=="middle"),])
fARMS <- row.names(metadata[(metadata$phase=="final"),])
baseARMS <- row.names(metadata[(metadata$phase=="initial") | (metadata$phase=="middle") | (metadata$phase=="final"),])

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

eutroARMS <- row.names(metadata[(metadata$treatment=="eutrophication"),])
fishARMS <- row.names(metadata[(metadata$treatment=="fishFarm"),])

wEutroARMS <- eutroARMS[c(eutroARMS %in% wARMS)]
eEutrpARMS <- eutroARMS[c(eutroARMS %in% eARMS)]

wFishARMS <- eutroARMS[c(fishARMS %in% wARMS)]
eFishARMS <- eutroARMS[c(fishARMS %in% eARMS)]

wEutroARMSrt <- wEutroARMS[c(wEutroARMS %in% RTARMS)]
eEutroARMSrt <- eEutrpARMS[c(eEutrpARMS %in% RTARMS)]

wEutroARMSrl <- wEutroARMS[c(wEutroARMS %in% RLARMS)]
eEutrpARMSrl <- eEutrpARMS[c(eEutrpARMS %in% RLARMS)]


# try to do a simply plot about eutro gain and lost
# 2nd year baseline 
wmBaseRich <- row.names(freqTable[rowSums((freqTable[,wmARMS]))>0,]) # 1765
emBaseRich <- row.names(freqTable[rowSums((freqTable[,emARMS]))>0,]) # 1802
mBaseRich <- union(wmBaseRich, emBaseRich) # 3009
mBaseRichShare <- wmBaseRich[wmBaseRich %in% emBaseRich] # 558

# eutro RT
wEutroRich <- row.names(freqTable[rowSums((freqTable[,wEutroARMSrt]))>0,]) # 1498
eEutroRich <- row.names(freqTable[rowSums((freqTable[,eEutroARMSrt]))>0,]) # 1612
EutroRich <- union(wEutroRich, eEutroRich) # 2290
eutroRichShare <- wEutroRich[wEutroRich %in% eEutroRich] # 820

# eutro lost anf eutro gain
wEutroLost <- wmBaseRich[!(wmBaseRich %in% wEutroRich)] # 1150/1765
wEutroGain <- wEutroRich[!(wEutroRich %in% wmBaseRich)] # 883/1498
wEutroShare <- wmBaseRich[(wmBaseRich %in% wEutroRich)] # 615

eEutroLost <- emBaseRich[!(emBaseRich %in% eEutroRich)] # 1100/1802
eEutroGain <- eEutroRich[!(eEutroRich %in% emBaseRich)] # 910/1612
eEutroShare <- emBaseRich[(emBaseRich %in% eEutroRich)] # 702


## okay so the above is th silly way now try make a table with
## for loop to get all the calculation out 

#################################

# resistance 
# build a plot data 
plotData <- as.data.frame(matrix(0, nrow=63, ncol=8))
names(plotData) <- c("phyla","treatment","water","baseLine","total","lost","gain","shared")
plotData$phyla <- c(rep("ALL",9),rep("Arthropoda",9),rep("Annelida",9),
                        rep("Bacillariophyta",9),rep("Rhodophyta",9),rep("Mollusca",9),rep("Porifera",9))


plotData$treatment <-rep(c("eutrophication","eutrophication","eutrophication",
                        "fishFarm","fishFarm","fishFarm",
                        "mBaseline","mBaseline","mBaseline"),7)


plotData$water <- rep(c("east", "west","combined"), 21)


freqTableALL <- freqTable
freqTableArthropoda <- freqTable[row.names(freqTable)%in%SOM[SOM$phylum == "Arthropoda",]$id,] 
freqTableAnnelida <- freqTable[row.names(freqTable)%in%SOM[SOM$phylum == "Annelida",]$id,] 
freqTableBacillariophyta <- freqTable[row.names(freqTable)%in%SOM[SOM$phylum == "Bacillariophyta",]$id,] 
freqTableRhodophyta <- freqTable[row.names(freqTable)%in%SOM[SOM$phylum == "Rhodophyta",]$id,] 
freqTableMollusca <- freqTable[row.names(freqTable)%in%SOM[SOM$phylum == "Mollusca",]$id,] 
freqTablePorifera <- freqTable[row.names(freqTable)%in%SOM[SOM$phylum == "Porifera",]$id,] 


### without mBaseline
for (k in c("ALL","Arthropoda","Annelida","Bacillariophyta","Rhodophyta","Mollusca","Porifera")) {
  freqTable <- get(paste0("freqTable",k))
  
  wiBaseRich <- row.names(freqTable[rowSums((freqTable[,wiARMS]))>0,]) # 1911
  eiBaseRich <- row.names(freqTable[rowSums((freqTable[,eiARMS]))>0,]) # 1867
  iBaseRich <- union(wiBaseRich, eiBaseRich) # 2849
  
  plotData[plotData$water=="west" & plotData$phyla == k,]$baseLine <- length(wiBaseRich)
  plotData[plotData$water=="east" & plotData$phyla == k,]$baseLine <- length(eiBaseRich)
  plotData[plotData$water=="combined" & plotData$phyla == k,]$baseLine <- length(iBaseRich)
  
  
  for (i in c("eutrophication","fishFarm")){
  for (j in c("west", "east")){
    wtargetRow <- row.names(metadata[(metadata$treatment==i & metadata$waterAct == j
                                      & metadata$water1=="west" & metadata$phase=="resistance"),])
    etargetRow <- row.names(metadata[(metadata$treatment==i & metadata$waterAct == j
                                      & metadata$water1=="east" & metadata$phase=="resistance"),])
    wbaseROW <- row.names(metadata[(metadata$treatment=="iBaseline" & metadata$water1=="west"),])
    ebaseROW <- row.names(metadata[(metadata$treatment=="iBaseline" & metadata$water1=="east"),])
    
    
    wspRich <- row.names(freqTable[(freqTable[,wtargetRow])>0,]) 
    espRich <- row.names(freqTable[(freqTable[,etargetRow])>0,]) 
    spRich <- union(wspRich,espRich)
    
    wbaseRich <- row.names(freqTable[rowSums((freqTable[,wbaseROW]))>0,])
    ebaseRich <- row.names(freqTable[rowSums((freqTable[,ebaseROW]))>0,])
    baseRich <- union(wbaseRich,ebaseRich)
      
    # the total 
    
    plotData[(plotData$treatment == i) & (plotData$water == j)& (plotData$phyla == k),]$total <-  length(spRich)
    
    # the lost 
    spLost <- baseRich[!(baseRich %in% spRich)] 
    plotData[(plotData$treatment == i) & (plotData$water == j)& (plotData$phyla == k),]$lost <-  length(spLost)
    
    # the gain 
    spGain <- spRich[!(spRich %in% baseRich)] 
    plotData[(plotData$treatment == i) & (plotData$water == j)& (plotData$phyla == k),]$gain <-  length(spGain)
    
    # the share 
    spShare <- spRich[(spRich %in% baseRich)] 
    plotData[(plotData$treatment == i) & (plotData$water == j)& (plotData$phyla == k),]$shared <-  length(spShare)
    
    
  }
  # the combined 
  targetRow1 <- row.names(metadata[(metadata$treatment==i & metadata$phase=="resistance"),])
  spRich1 <- row.names(freqTable[rowSums((freqTable[,targetRow1]))>0,]) 
  plotData[(plotData$treatment == i) & (plotData$water == "combined") 
           & (plotData$phyla == k),]$total <-  length(spRich1)

  # the combined lost 
  baseRich1 <- union(eiBaseRich,wiBaseRich)
  spLost1 <- baseRich1[!(baseRich1 %in% spRich1)] 
  plotData[(plotData$treatment == i) & (plotData$water == "combined") & (plotData$phyla == k),]$lost <-  length(spLost1)

  # the gain 
  spGain1 <- spRich1[!(spRich1 %in% baseRich1)] 
  plotData[(plotData$treatment == i) & (plotData$water == "combined") & (plotData$phyla == k),]$gain <-  length(spGain1)
  
  # the share 
  spShare1 <- spRich1[(spRich1 %in% baseRich1)] 
  plotData[(plotData$treatment == i) & (plotData$water == "combined") & (plotData$phyla == k),]$shared <-  length(spShare1)
  
  }

}

### following only for mBaseline
for (k in c("ALL","Arthropoda","Annelida","Bacillariophyta","Rhodophyta","Mollusca","Porifera")) {
  freqTable <- get(paste0("freqTable",k))
  
  wiBaseRich <- row.names(freqTable[rowSums((freqTable[,wiARMS]))>0,]) # 1911
  eiBaseRich <- row.names(freqTable[rowSums((freqTable[,eiARMS]))>0,]) # 1867
  iBaseRich <- union(wiBaseRich, eiBaseRich) # 2849


      wtargetRow <- row.names(metadata[(metadata$treatment== "mBaseline" & metadata$water1=="west"),])
      etargetRow <- row.names(metadata[(metadata$treatment== "mBaseline" & metadata$water1=="east"),])
      wbaseROW <- row.names(metadata[(metadata$treatment=="iBaseline" & metadata$water1=="west"),])
      ebaseROW <- row.names(metadata[(metadata$treatment=="iBaseline" & metadata$water1=="east"),])
      
      
      wspRich <- row.names(freqTable[(freqTable[,wtargetRow])>0,]) 
      espRich <- row.names(freqTable[(freqTable[,etargetRow])>0,]) 
      spRich <- union(wspRich,espRich)
      
      wbaseRich <- row.names(freqTable[rowSums((freqTable[,wbaseROW]))>0,])
      ebaseRich <- row.names(freqTable[rowSums((freqTable[,ebaseROW]))>0,])
      baseRich <- union(wbaseRich,ebaseRich)
      
      # the total 
      
      plotData[(plotData$treatment == "mBaseline") & (plotData$water == "combined")& (plotData$phyla == k),]$total <-  length(spRich)
      
      # the lost 
      spLost <- baseRich[!(baseRich %in% spRich)] 
      plotData[(plotData$treatment == "mBaseline") & (plotData$water == "combined")& (plotData$phyla == k),]$lost <-  length(spLost)
      
      # the gain 
      spGain <- spRich[!(spRich %in% baseRich)] 
      plotData[(plotData$treatment == "mBaseline") & (plotData$water == "combined")& (plotData$phyla == k),]$gain <-  length(spGain)
      
      # the share 
      spShare <- spRich[(spRich %in% baseRich)] 
      plotData[(plotData$treatment == "mBaseline") & (plotData$water == "combined")& (plotData$phyla == k),]$shared <-  length(spShare)
      
      
    }
  

# only look at the combined data because they are the only on that matter
plotData.fi <- plotData %>% filter(water =="combined")


## play with plotData
plotData.fi$GvL <- plotData.fi$gain - plotData.fi$lost
plotData.fi$GvLper <- 100 * (plotData.fi$GvL / plotData.fi$baseLine)
plotData.fi$deeds <- paste0(plotData.fi$phyla,"_", plotData.fi$treatment,"_",plotData.fi$water)
plotData.fi$deeds1 <- paste0(plotData.fi$water,"_",plotData.fi$phyla)
plotData.fi$deeds2 <- paste0(plotData.fi$treatment,"_",plotData.fi$water )
write.csv(plotData.fi,"decontamByArms/noSED/sidewaybar_2years.csv")

# do some stats 
model1 <- aov(GvL ~ treatment, data = plotData.fi%>% filter((phyla != "ALL") & (water == "combined")))
TukeyHSD(model1) # eutrophication significantlly lose 

model2 <- aov(GvL ~ phyla+treatment, data = plotData.fi%>% filter((phyla != "ALL") & (water == "combined")))
TukeyHSD(model2)


# make some plots 
# make long list to plot
data_long <- pivot_longer(plotData.fi, cols = c(lost, gain), names_to = "gain_lost", 
                          values_to = "OTUs")

# negative the lost value 
data_long[data_long$`gain_lost`=="lost",]$OTUs <- 
  data_long[data_long$`gain_lost`=="lost",]$OTUs * -1

# percentage the lost/gain 
data_long$perGL <- round((data_long$GvL/data_long$baseLine) * 100, 1)
data_long$perGLS <- -23
data_long[data_long$gain_lost=="gain",]$perGLS <- round((data_long%>%filter(gain_lost=="gain") %>% .$OTUs/
                                                          data_long%>%filter(gain_lost=="gain") %>% .$baseLine) * 100, 1)
data_long[data_long$gain_lost=="lost",]$perGLS <- round((data_long%>%filter(gain_lost=="lost") %>% .$OTUs/
                                                          data_long%>%filter(gain_lost=="lost") %>% .$baseLine) * 100, 1)

# Create the horizontal bar plot


sidebarALLbyALL_GLS <- ggplot(data_long%>%filter(phyla =="ALL"& water == "combined"), aes(x = perGLS, y = deeds1, fill = gain_lost)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +
  scale_fill_manual(values  = c("gain" = "#0072B2", "lost" = "#D55E00")) +
  geom_text(data = subset(data_long%>%filter(phyla =="ALL"& water == "combined"), gain_lost == "gain"), aes(label = perGLS), 
            position = position_stack(vjust = 0), hjust=-0.5, color = "yellow", size = 5) +  # Center "gain" labels at the top
  geom_text(data = subset(data_long%>%filter(phyla =="ALL"& water == "combined"), gain_lost == "lost"), aes(label = perGLS), 
            position = position_stack(vjust = 0), hjust=-0.5,color = "black", size = 5) +  # Center "gain" labels at the top
  
  #labs(title = "Horizontal Bar Plot with Different Treatments",x = "OTUs", y = "Treatment") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none")+
  labs(x = NULL, y = NULL) +
  facet_wrap(~treatment)

data_long$phyla <- factor(data_long$phyla, levels= c("ALL", "Arthropoda","Annelida",
                          "Bacillariophyta", "Rhodophyta", "Mollusca", "Porifera")) # reoder the y axis 
data_long$water <- factor(data_long$water, levels= c("west", "east","combined")) # reoder the y axis 


data_long$phyla <- factor(data_long$phyla, levels = c("Porifera","Rhodophyta","Mollusca","Bacillariophyta","Annelida","Arthropoda","ALL"))
sidebarT6byTreatment_GLS <- ggplot(data_long%>%filter(phyla !="ALL" & water == "combined"), aes(x = perGLS, y = phyla, fill = gain_lost)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +
  scale_fill_manual(values  = c("gain" = "#0072B2", "lost" = "#D55E00")) +
  geom_text(data = subset(data_long%>%filter(phyla !="ALL"& water == "combined"), gain_lost == "gain"), aes(label = perGLS), 
            position = position_stack(vjust = 0), hjust=-0.75, color = "yellow", size = 3.5) +  # Center "gain" labels at the top
  geom_text(data = subset(data_long%>%filter(phyla !="ALL"& water == "combined"), gain_lost == "lost"), aes(label = perGLS), 
            position = position_stack(vjust = 0), hjust=-0.1,color = "black", size = 3.5) +  # Center "gain" labels at the top
  #labs(title = "Horizontal Bar Plot with Different Treatments",x = "OTUs", y = "Treatment") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none")+
  labs(x = NULL, y = NULL) +
  facet_wrap(~ interaction(treatment)) 


data_long$treatment <- factor(data_long$treatment, levels= c("eutrophication","fishFarm", "mBaseline")) # reoder the treatment 

data_long$phyla <- factor(data_long$phyla, levels = c("Arthropoda","Annelida","Bacillariophyta","Mollusca","Rhodophyta","Porifera","ALL"))
allLinePlot <- ggplot(data_long%>%filter(water == "combined" & gain_lost=="gain"), aes(x=treatment, y=GvL/baseLine, group=phyla)) +
  geom_point(aes(color = phyla, shape = phyla), size = 3,
             position = position_jitter(width = 0.2, height = 0)  # Horizontal jitter only
  ) + 
  scale_shape_manual(values = c(16, 17, 18, 15, 3, 4, 8)) +
theme_minimal() +
  ylab("")+
  xlab("")+
  scale_x_discrete(labels = c(
    "eutrophication" = "Domestic Sewage",
    "fishFarm" = "Mariculture",
    "mBaseline" = "Marine Park"
  )) +
  theme(
    legend.title = element_blank()  # This removes the legend title
  )



combined_plot <- plot_grid(sidebarALLbyALL_GLS,sidebarT6byTreatment_GLS,allLinePlot, ncol = 1) # plot the above together 


####################################
