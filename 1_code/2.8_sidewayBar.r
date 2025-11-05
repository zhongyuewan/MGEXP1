#########################
# PREPARE R ENVIRONMENT 
#########################
library(Biostrings)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(forcats)
library(tidyr)
library(cowplot)
library(scales)



## set WD 
setwd("/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#1 biodiveristy/#2 CRF EXP1/#2 dataAnalysis_all/1_data")
rm(list=ls())

# get some color see see
myColor <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#79bfe7","#D55E00","#CC79A7")

## get data
freqTable <- read.csv('decontamByArms/noSED/freqTable.csv', header = TRUE, row.names = 1)
metadata <- read.csv('decontamByArms/noSED/metadata.csv', row.names = 1)
sequenceTable <- readDNAStringSet('decontamByArms/noSED/dna-sequences97.fasta')
SOM.d1 <- read.csv('taxAssign/blast/TaxAsn_shelbyOmidori.csv') 

SOM <- SOM.d1 %>% filter(id %in% row.names(freqTable))
nrow(SOM.d1) - nrow(SOM) # 408 sedimentation rows 

plotData.fi.rt <- read.csv("decontamByArms/noSED/sidewaybar_2years.csv")
plotData.fi.rl <- read.csv("decontamByArms/noSED/sidewaybar_3years.csv")

# fit formate for plot 
plotData.fi.rt[1:3,]$deeds1 <- "Total Richness"
plotData.fi.rl[1:3,]$deeds1 <- "Total Richness"


plotData.fi.rt$treatment <- factor(plotData.fi.rt$treatment, levels = c("mBaseline","fishFarm","eutrophication"))
plotData.fi.rl$treatment <- factor(plotData.fi.rl$treatment, levels = c("mBaseline","fishFarm","eutrophication"))


# make some plots 
# make long list to plot
data_long.rt <- pivot_longer(plotData.fi.rt, cols = c(lost, gain), names_to = "gain_lost", 
                          values_to = "OTUs")

data_long.rl <- pivot_longer(plotData.fi.rl, cols = c(lost, gain), names_to = "gain_lost", 
                             values_to = "OTUs")


# negative the lost value 
data_long.rt[data_long.rt$`gain_lost`=="lost",]$OTUs <- 
  data_long.rt[data_long.rt$`gain_lost`=="lost",]$OTUs * -1

data_long.rl[data_long.rl$`gain_lost`=="lost",]$OTUs <- 
  data_long.rl[data_long.rl$`gain_lost`=="lost",]$OTUs * -1



# percentage the lost/gain 
data_long.rt$perGL <- round((data_long.rt$GvL/data_long.rt$baseLine) * 100, 1)
data_long.rl$perGL <- round((data_long.rl$GvL/data_long.rl$baseLine) * 100, 1)


data_long.rt$perGLS <- -23 # just assign a random number to start a column 
data_long.rt[data_long.rt$gain_lost=="gain",]$perGLS <- round((data_long.rt%>%filter(gain_lost=="gain") %>% .$OTUs/
                                                          data_long.rt%>%filter(gain_lost=="gain") %>% .$baseLine) * 100, 1)
data_long.rt[data_long.rt$gain_lost=="lost",]$perGLS <- round((data_long.rt%>%filter(gain_lost=="lost") %>% .$OTUs/
                                                          data_long.rt%>%filter(gain_lost=="lost") %>% .$baseLine) * 100, 1)


data_long.rl$perGLS <- -23 # just assign a random number to start a column 
data_long.rl[data_long.rl$gain_lost=="gain",]$perGLS <- round((data_long.rl%>%filter(gain_lost=="gain") %>% .$OTUs/
                                                                 data_long.rl%>%filter(gain_lost=="gain") %>% .$baseLine) * 100, 1)
data_long.rl[data_long.rl$gain_lost=="lost",]$perGLS <- round((data_long.rl%>%filter(gain_lost=="lost") %>% .$OTUs/
                                                                 data_long.rl%>%filter(gain_lost=="lost") %>% .$baseLine) * 100, 1)

# Create the horizontal bar plot
# for RT

sidebarALLbyALL_GLS.rt <- ggplot(data_long.rt%>%filter(phyla =="ALL"), aes(x = perGLS, y = deeds1, fill = gain_lost)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +
  scale_fill_manual(values  = c("gain" = "#79bfe7", "lost" = "#D55E00")) +
  scale_x_continuous(limits = c(-70, 70)) +
  geom_text(data = subset(data_long.rt%>%filter(phyla =="ALL"), gain_lost == "gain"), aes(label = number(perGLS, accuracy = 0.1)), 
            position = position_stack(vjust=0.75), hjust=1, color = "black", size = 4) +  
  geom_text(data = subset(data_long.rt%>%filter(phyla =="ALL"), gain_lost == "lost"), aes(label = number(perGLS, accuracy = 0.1)), 
            position = position_stack(vjust=0.25), hjust=0, color = "black", size = 4) +  
  
  #labs(title = "Horizontal Bar Plot with Different Treatments",x = "OTUs", y = "Treatment") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 13),)+
  labs(x = NULL, y = NULL) +
  facet_wrap(~treatment)+
  theme(
    strip.text = element_blank()
  )

data_long.rt$phyla <- factor(data_long.rt$phyla, levels= c("ALL","ALLIDed", "Arthropoda","Annelida",
                          "Bacillariophyta", "Rhodophyta", "Mollusca", "Porifera")) # reoder the y axis 



data_long.rt$phyla <- factor(data_long.rt$phyla, levels = c("Porifera","Rhodophyta","Mollusca","Bacillariophyta","Annelida","Arthropoda","ALLIDed","ALL"))
sidebarT6byTreatment_GLS.rt <- ggplot(data_long.rt%>%filter(phyla !="ALL" & phyla != "ALLIDed"), aes(x = perGLS, y = phyla, fill = gain_lost)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +
  scale_fill_manual(values  = c("gain" = "#79bfe7", "lost" = "#D55E00")) +
  scale_x_continuous(limits = c(-150, 155)) +
  geom_text(data = subset(data_long.rt%>%filter(phyla !="ALL" & phyla != "ALLIDed"), gain_lost == "gain"), aes(label = number(perGLS, accuracy = 0.1)), 
            position = position_stack(vjust=1), hjust=0, color = "black", size = 4) +  
  geom_text(data = subset(data_long.rt%>%filter(phyla !="ALL" & phyla != "ALLIDed"), gain_lost == "lost"), aes(label = number(perGLS, accuracy = 0.1)), 
            position = position_stack(vjust=0), hjust=1, color = "black", size = 4) +  
  #labs(title = "Horizontal Bar Plot with Different Treatments",x = "OTUs", y = "Treatment") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 13),)+
  labs(x = NULL, y = NULL) +
  facet_wrap(~ interaction(treatment)) +
  theme(
    strip.text = element_blank()
  )



# for RL
sidebarALLbyALL_GLS.rl <- ggplot(data_long.rl%>%filter(phyla =="ALL"), aes(x = perGLS, y = deeds1, fill = gain_lost)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +
  scale_fill_manual(values  = c("gain" = "#79bfe7", "lost" = "#D55E00")) +
  scale_x_continuous(limits = c(-70,70)) +
  geom_text(data = subset(data_long.rl%>%filter(phyla =="ALL"), gain_lost == "gain"), aes(label = number(perGLS, accuracy = 0.1)), 
            position = position_stack(vjust=0.75), hjust=1, color = "black", size = 4) +  
  geom_text(data = subset(data_long.rl%>%filter(phyla =="ALL"), gain_lost == "lost"), aes(label = number(perGLS, accuracy = 0.1)), 
            position = position_stack(vjust=0.25), hjust=0, color = "black", size = 4) +  
  
  #labs(title = "Horizontal Bar Plot with Different Treatments",x = "OTUs", y = "Treatment") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 13),)+
  labs(x = NULL, y = NULL) +
  facet_wrap(~treatment)+
  theme(
    strip.text = element_blank()
  )

data_long.rl$phyla <- factor(data_long.rl$phyla, levels= c("ALL","ALLIDed", "Arthropoda","Annelida",
                                                           "Bacillariophyta", "Rhodophyta", "Mollusca", "Porifera")) # reoder the y axis 



data_long.rl$phyla <- factor(data_long.rl$phyla, levels = c("Porifera","Rhodophyta","Mollusca","Bacillariophyta","Annelida","Arthropoda","ALLIDed","ALL"))
sidebarT6byTreatment_GLS.rl <- ggplot(data_long.rl%>%filter(phyla !="ALL" & phyla != "ALLIDed"), aes(x = perGLS, y = phyla, fill = gain_lost)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +
  scale_fill_manual(values  = c("gain" = "#79bfe7", "lost" = "#D55E00")) +
  scale_x_continuous(limits = c(-150, 155)) +
  geom_text(data = subset(data_long.rl%>%filter(phyla !="ALL" & phyla != "ALLIDed"), gain_lost == "gain"), 
            aes(label = number(perGLS, accuracy = 0.1)), 
            position = position_stack(vjust=1), hjust=0, color = "black", size = 4) +  
  geom_text(data = subset(data_long.rl%>%filter(phyla !="ALL" & phyla != "ALLIDed"), gain_lost == "lost"), 
            aes(label = number(perGLS, accuracy = 0.1)), 
            position = position_stack(vjust=0), hjust=1,color = "black", size = 4) + 
  #labs(title = "Horizontal Bar Plot with Different Treatments",x = "OTUs", y = "Treatment") +
  theme_minimal() +
  theme(axis.line = element_line(color = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",    
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 13),)+
  labs(x = NULL, y = NULL) +
  facet_wrap(~ interaction(treatment)) +
  theme(
    strip.text = element_blank()
  )




combined_plot <- plot_grid(sidebarALLbyALL_GLS.rt,sidebarT6byTreatment_GLS.rt,sidebarALLbyALL_GLS.rl,sidebarT6byTreatment_GLS.rl, ncol = 1) # plot the above together 




#### following gonna do chi square 
data.chi.gain.rt.d1 <- (as.data.frame(data_long.rt) %>%  filter(gain_lost=="gain"))[,c(2,3,14)]
data.chi.gain.rt <- as.data.frame(data.chi.gain.rt.d1 %>%
  pivot_wider(
    names_from = treatment,    # Column to get new column names from
    values_from = OTUs         # Column to get values from
  ))

data.chi.lost.rt.d1 <- (as.data.frame(data_long.rt) %>%  filter(gain_lost=="lost"))[,c(2,3,14)]
data.chi.lost.rt <- as.data.frame(data.chi.lost.rt.d1 %>%
  pivot_wider(
    names_from = treatment,    # Column to get new column names from
    values_from = OTUs         # Column to get values from
  ))

model1 <- chisq.test(data.chi.gain.rt[3:8,])
model1 # sig
model1$stdres

model2 <- chisq.test(-data.chi.lost.rt[3:8,-1])
model2 # not sig 
model2$stdres

data.chi.net.rt <- data.chi.gain.rt + data.chi.lost.rt
data.chi.net.rt[,1] <- data.chi.gain.rt[,1]


data.chi.gain.rl.d1 <- (as.data.frame(data_long.rl) %>%  filter(gain_lost=="gain"))[,c(2,3,13)]
data.chi.gain.rl <- as.data.frame(data.chi.gain.rl.d1 %>%
                                    pivot_wider(
                                      names_from = treatment,    # Column to get new column names from
                                      values_from = OTUs         # Column to get values from
                                    ))

data.chi.lost.rl.d1 <- (as.data.frame(data_long.rl) %>%  filter(gain_lost=="lost"))[,c(2,3,13)]
data.chi.lost.rl <- as.data.frame(data.chi.lost.rl.d1 %>%
                                    pivot_wider(
                                      names_from = treatment,    # Column to get new column names from
                                      values_from = OTUs         # Column to get values from
                                    ))

model3 <- chisq.test(data.chi.gain.rl[3:8,-1])
model3 # sig
model3$stdres

model4 <- chisq.test(-data.chi.lost.rl[3:8,-1])
model4 # sig 
model4$stdres



