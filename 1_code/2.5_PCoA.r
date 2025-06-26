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

## get data
# it could be jaccard or BC.H, just load the differernt diversity matrix 
# BC.H here is not really BC h, it is BC on square root transfered count data
freqTable.sed <- read.csv('decontamByArms/freqTableCleanByArms.csv', header = TRUE, row.names = 1)
metadata.sed <- read.csv('decontamByArms/sample-metadata.byArms.csv', row.names = 1)

## fix table names and stuffs 
names(freqTable.sed) <- gsub("X","", names(freqTable.sed))

nosedARMS <- row.names(metadata.sed[metadata.sed$treatment!="sedimentation",])
metadata.sed$arms <- row.names(metadata.sed)
metadata <- metadata.sed %>% filter (arms %in% nosedARMS)
row.names(metadata) <- metadata$arms
freqTable.d1 <- freqTable.sed[,nosedARMS] 
freqTable <- freqTable.d1 %>% filter (rowSums(.)!= 0) # 7696 

# make distance matrix
distMeuc <- as.matrix(vegdist(t(freqTable), method = "euclidean")) 
distMBC <-as.matrix(vegdist(t(freqTable), method = "bray")) 
distMBCsr <-as.matrix(vegdist(t(sqrt(freqTable)), method = "bray")) 
distMJC <-as.matrix(vegdist(t(freqTable), method = "jaccard")) 

metadata$time <- "2nd Year" # add a new column for time 
metadata[metadata$dayTotal<400,]$time <- "1st Year" # add a new column for time 
metadata[metadata$dayTotal>800,]$time <- "3rd Year" # add a new column for time 
metadata[metadata$treatment=="iBaseline" | 
           metadata$treatment=="mBaseline" |
           metadata$treatment=="fBaseline",]$treatment <- "Baseline" # add a new column for time 

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


###### plot succession 
metadata$ARMS <- row.names(metadata)  
plot5ARMS <- metadata%>% filter(treatment=="Baseline")%>% .$ARMS
dist_matrix <- prcomp(distMBCsr[plot5ARMS,plot5ARMS]) 

  evalue <- fviz_eig(dist_matrix, addlabels = TRUE) 
  dataPlot <- data.frame(matrix(nrow=length(plot5ARMS), ncol=4))
  row.names(dataPlot) <- plot5ARMS
  colnames(dataPlot) <- c("PC1", "PC2", "treatment", "water")
  
  dataPlot$PC1 <- dist_matrix$x[,1]
  dataPlot$PC2 <- dist_matrix$x[,2]
  dataPlot$treatment <- metadata[plot5ARMS,]$treatment
  dataPlot$water1 <- metadata[plot5ARMS,]$water1
  dataPlot$waterAct <- metadata[plot5ARMS,]$waterAct
  dataPlot$time <- metadata[plot5ARMS,]$time
  dataPlot$timewaterAct <- paste0(metadata[plot5ARMS,]$time, metadata[plot5ARMS,]$waterAct)
  dataPlot$timewater1 <- paste0(metadata[plot5ARMS,]$time, metadata[plot5ARMS,]$water1)

  plotSuccession  <- ggplot(dataPlot, aes(x = PC1, y = PC2, color = water1, shape = time)) +
    geom_point(aes(shape = time,size = 80)) +
    scale_color_manual(values = c(
      "east" = "royalblue",   
      "west" = "indianred1")) +
    xlab("PC1") +
    ylab("PC2") +
    theme(text = element_text(size = 35)) +
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank())+
    
    stat_ellipse(
      aes(group = time),  # Group ellipses by `water1`
      level = 0.95,         # 95% confidence level (adjust as needed)
      linewidth = 0.5,        # Line thickness
      linetype = 1   # Line style
    ) +
    
    theme_classic()
  
  
  
###### plot all
  plot4ARMS_nosed <- metadata%>% filter(treatment!="sedimentation")%>% .$ARMS
  
  dist_matrix <- prcomp(distMBCsr[plot4ARMS_nosed,plot4ARMS_nosed]) # bray crutit with squart root which is the best to visualization 
  # dist_matrix <- prcomp(distMJC[plot4ARMS,plot4ARMS]) #jaccard dist
  # dist_matrix <- prcomp(distMeuc[plot4ARMS,plot4ARMS]) # euclidean dist 
  # dist_matrix <- prcomp(distMBC[plot4ARMS,plot4ARMS])  # bray crutit 
  
  
  evalue <- fviz_eig(dist_matrix, addlabels = TRUE) # 34.4%, 12.1%
  dataPlot <- data.frame(matrix(nrow=length(plot4ARMS_nosed), ncol=4))
  row.names(dataPlot) <- plot4ARMS_nosed
  colnames(dataPlot) <- c("PC1", "PC2", "treatment", "water")
  
  dataPlot$PC1 <- dist_matrix$x[,1]
  dataPlot$PC2 <- dist_matrix$x[,2]
  dataPlot$treatment <- metadata[plot4ARMS_nosed,]$treatment
  dataPlot$water1 <- metadata[plot4ARMS_nosed,]$water1
  dataPlot$waterAct <- metadata[plot4ARMS_nosed,]$waterAct
  dataPlot$time <- metadata[plot4ARMS_nosed,]$time
  dataPlot$timewaterAct <- paste0(metadata[plot4ARMS_nosed,]$time, metadata[plot4ARMS_nosed,]$waterAct)
  dataPlot$timewater1 <- paste0(metadata[plot4ARMS_nosed,]$time, metadata[plot4ARMS_nosed,]$water1)
  
  
  dataPlot$timewaterAct <- factor(dataPlot$timewaterAct, levels = c("1st Yeareast", "2nd Yeareast", "3rd Yeareast",
         "1st Yearwest", "2nd Yearwest", "3rd Yearwest"))
  dataPlot$timewater1 <- factor(dataPlot$timewater1, levels = c("1st Yeareast", "2nd Yeareast", "3rd Yeareast",
                                                                    "1st Yearwest", "2nd Yearwest", "3rd Yearwest"))

  
# try to do it with shleby's idea 
  # Create the revised plot with solid and hollow symbols
plotALL  <- ggplot(dataPlot, aes(x = PC1, y = PC2, color = treatment, shape = timewaterAct)) +
    geom_point(aes(shape = timewaterAct,
                   size = 80,
                   stroke=2)) +
    scale_shape_manual(values = c(1,2,0,16, 17, 15))+
    scale_color_manual(values = c(
    "Baseline" = "#0072B2",   
    "fishFarm" = "#E69F00",   
    "eutrophication" = "#D55E00")) +
    xlab("PC1 (34.4%)") +
    ylab("PC2 (12.1%)") +
    theme_classic()+
    theme(text = element_text(size = 17)) +
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          legend.position = "none")



