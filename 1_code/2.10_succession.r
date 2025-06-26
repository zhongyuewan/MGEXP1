#########################
# PREPARE R ENVIRONMENT #
#########################
library(Biostrings)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(forcats)
library(cowplot)
library(MASS)

## set WD 
setwd("/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#1 biodiveristy/#2 CRF EXP1/#2 dataAnalysis_all/1_data")
rm(list=ls())

## get data
freqTable <- read.csv('decontamByArms/noSED/freqTable.csv', header = TRUE, row.names = 1)
metaData <- read.csv('decontamByArms/noSED/metadata.csv', row.names = 1)
sequenceTable <- readDNAStringSet('decontamByArms/noSED/dna-sequences97.fasta')
SOM.d1 <- read.csv('taxAssign/blast/TaxAsn_shelbyOmidori.csv') 


SOM <- SOM.d1 %>% filter(id %in% row.names(freqTable))
nrow(SOM.d1) - nrow(SOM) # 408 sedimentation rows 



# summary of the above
# shelby DB is able to pick up 3043 sequnce passing the filter 
# but midor has been able to pick up more with higher score so it out perfor shelby db
# but because shelby is a locally curated DB so i will use shelby over midor "SOM"
# which are the shared pick up (3444) + shared picked up bias towards shelby (973) = 4417 

#add a $year to show 1year, 2years, 3years arms 
metaData$years <- "Resilience"
metaData[metaData$phase=="initial",]$years <- "Seeding"
metaData[metaData$phase=="middle" | metaData$phase=="resistance" ,]$years <- "Resistance"
metaData$years <- factor(metaData$years, levels = c("Seeding", "Resistance", "Resilience"))

#touch up table 

names(freqTable) <- gsub("X","", names(freqTable))


freqTable$phylum <- "unAssigned" 
##### the following code is horriblly wrong but i am happy i found it out ##########
##### no i need to redo tax assign, hopfully the data is not off by far ############
## freqTable[rownames(freqTable) %in% SOM$id,]$phylum <- SOM$phylum #dumbbbbb ######
################# the above is a dumb code that i thought at the time was smart ####
## it does NOT align the ids with the phylum so it just copy/paste the phylum -->###
## to the freqTable with a random order so the seq is no longer matching the phy ### 

# so i need to build a forloop to make sure the id names are aligned 
for (i in 1: nrow(SOM)) {
  SeqName <- SOM$id[i]
  SeqPhy <- SOM$phylum[i]
  
  freqTable[rownames(freqTable)==SeqName,]$phylum <- SeqPhy
  
}




rankPhylum <- function (df){
  
  result <- df %>%
    group_by(phylum) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    mutate(rank = row_number())
  
  return(result)
  
}


# make a new long df to store the data for a plot 
temp <- freqTable[freqTable[,1] >0,][,c(1,29)] # the first arms 
temp <- rankPhylum(temp) 

plot_data <- as.data.frame(matrix(0,nrow=nrow(temp), ncol=3))
names(plot_data) <- c("ARMS","phylum","count")                           
           
# wirte the first in and the loop it                 
plot_data$ARMS <- names(freqTable)[1]
temp <- freqTable[freqTable[,1] >0,][,c(1,29)] # the first arms 
plot_data$phylum <- rankPhylum(temp)$phylum
plot_data$count <- rankPhylum(temp)$count


for (i in 2:28){
  
  # count the rank
  temp <- freqTable[freqTable[,i] >0,][,c(i,29)] # the ith arms 
  new_data <- as.data.frame(matrix(0, nrow = nrow(rankPhylum(temp)), ncol=3))
  names(new_data) <- c("ARMS","phylum","count")   
  
  # write data in loop
  new_data$ARMS <- names(freqTable)[i]
  new_data$phylum <- rankPhylum(temp)$phylum
  new_data$count <- rankPhylum(temp)$count
  
  # loop it 
  
  plot_data <- rbind(plot_data, new_data)
  
  
}

#make a function to add a all phylum count 
addCount <- function (df){
  temp <- as.data.frame(matrix(0, nrow=28, ncol=3))
  names(temp) <- c("ARMS","phylum","count")
  row.names(temp) <- names(freqTable)[1:28]
  
  for (i in 1:28){
  temp[i,1] <- names(freqTable)[i]
  temp[i,2] <- "allcount"
  temp[i,3] <- sum(plot_data$count[plot_data$ARMS==names(freqTable[i])])
}
  return(temp)
  
}

all_data <- addCount(plot_data)
plot_data <- rbind(plot_data,all_data)


# rearrange the plot for better reference 
# add treatment/phase/water 
getColumn <- function(df,extract){
  for (i in 1:nrow(df)) {
    ARMSn <- df[i,"ARMS"]
    df[i,extract] <- metaData[ARMSn, extract]
    
  }
  return(df)
}

# add in more colums from metadata 
plot_data1 <- getColumn(plot_data,"waterAct")
plot_data1 <- getColumn(plot_data1,"treatment")
plot_data1 <- getColumn(plot_data1,"phase")
plot_data1 <- getColumn(plot_data1,"dayTotal")
plot_data1 <- getColumn(plot_data1,"years")




# let's plot 
# find a good color
display.brewer.all() # look at the color options 
phylum_colors <- rep(brewer.pal(9, "Set1"), 6)# get nice color, set1 is nice 

# reorder y axis to rank from top to bottom 
plot_data1$phylum_ordered <- reorder(plot_data$phylum, -plot_data$count) # reorder data 

# reorder x axis to rank from initial, mid, final, RT, RL and the west to east. 
plot_data1$treatment <- factor(plot_data1$treatment, levels = c("iBaseline","mBaseline",
                                                                "fBaseline","sedimentation","fishFarm", "eutrophication"))
plot_data1$waterAct <- factor(plot_data1$waterAct, levels = c("west","east"))
plot_data1$phase <- factor(plot_data1$phase, levels = c("initial","middle","resistance", "final","resilience"))

# what to rank? 
plot_data1$ARMS_reorder <- reorder(plot_data1$ARMS, as.numeric(plot_data1$waterAct)) # west/east
plot_data1$ARMS_reorder <- reorder(plot_data1$ARMS_reorder, as.numeric(plot_data1$treatment)) # treatment
plot_data1$ARMS_reorder <- reorder(plot_data1$ARMS_reorder, as.numeric(plot_data1$phase)) # phase


# plot 
plot1 <- ggplot(subset(plot_data1,phylum!="unAssigned" & phylum!="allcount"), aes(fill=phylum_ordered, y=count, x=ARMS_reorder)) + # bar plot by counts SKIP unAssigned
  geom_bar(position="stack", stat="identity")+
  guides(fill = guide_legend(order = 2))+
  scale_fill_manual(values = phylum_colors)

# move the "unAssigned" till the end
plot_data1$phylum_ordered <- fct_relevel(plot_data1$phylum_ordered, "unAssigned", after = Inf) # it worked 
plot2 <- ggplot(subset(plot_data1,phylum!="allcount"), aes(fill=phylum_ordered, y=count, x=ARMS_reorder)) + # bar plot by percentage 
  geom_bar(position="fill", stat="identity")+
  guides(fill = guide_legend(order = 1))+
  scale_fill_manual(values = phylum_colors)




plot.artho <- ggplot(plot_data1 %>% filter(phylum=="Arthropoda"),
                     aes(y=count, x=years, color=waterAct))+
  geom_boxplot()

plot.anne <- ggplot(plot_data1 %>% filter(phylum=="Annelida"),
                    aes(y=count, x=years, color=waterAct))+
  geom_boxplot()

plot.mollu <- ggplot(plot_data1 %>% filter(phylum=="Mollusca") ,
                     aes(y=count, x=years, color=waterAct))+
  geom_boxplot()


plot.Base <- ggplot(plot_data1 %>%filter(grepl("Baseline", treatment))%>% filter(phylum=="allcount") ,
                   aes(y=count, x=years, color=waterAct))+
  geom_boxplot()+
  theme_classic()


plot.all <- ggplot(plot_data1 %>% filter(phylum=="allcount") ,
                     aes(y=count, x=years, color=waterAct))+
  geom_boxplot()+
  theme_classic()




checkLMmodel <- function (df,p){
  
  model <- aov(count~years, df %>% filter(phylum==p))
  summary(model)
}
checkLMmodel(plot_data1, "Arthropoda") # p < 0.001 
checkLMmodel(plot_data1, "Annelida") # p < 0.001
checkLMmodel(plot_data1, "Mollusca") # p < 0.001

# turn phase into 12/24/30 months 
plot_data1$soaktime <- 12
plot_data1[plot_data1$phase=="resistance" | plot_data1$phase=="middle", ]$soaktime <- 24
plot_data1[plot_data1$phase=="resilience" | plot_data1$phase=="final", ]$soaktime <- 30


## gonna plot the baseline arms and see how they go 
## make treatment 
# this is the normal group 
plot_data1$waterAct <- factor(plot_data1$waterAct, levels = c("west","east"))
# make a forloop to make plots and get the stats 
# first creat a df to store the stats output 

succOutput <- as.data.frame(matrix(0, nrow=14,ncol=5))
succOutput[,1] <- rep(c("allcount", "Arthropoda", "Annelida", "Mollusca", "Rhodophyta","Porifera", "Bacillariophyta"),2)
names(succOutput) <- c("phylum","water","slope","zslope","pslope")
succOutput[1:7,2] <- "west"
succOutput[8:14,2] <- "east"

# 1 all 
for (i in c("allcount", "Arthropoda", "Annelida", "Mollusca", "Rhodophyta","Porifera", "Bacillariophyta")){
model_east <- glm.nb(count ~ soaktime, # + I(soaktime^2), 
                 # family = poisson,
                  data = plot_data1 %>% 
                    filter(phylum == i,
                           treatment %in% c("iBaseline", "mBaseline", "fBaseline"),
                           waterAct == "east"))

model_west <- glm.nb(count ~ soaktime, # + I(soaktime^2), 
                # family = poisson,
                  data = plot_data1 %>% 
                    filter(phylum == i,
                           treatment %in% c("iBaseline", "mBaseline", "fBaseline"),
                           waterAct == "west"))
# east 
pred_east <- data.frame(
  soaktime = seq(min(plot_data1$soaktime), 
                 max(plot_data1$soaktime), 
                 length.out = 100),
  waterAct = "east"
)
pred_east$count <- predict(model_east, newdata = pred_east, type = "response")

# west 
pred_west <- data.frame(
  soaktime = seq(min(plot_data1$soaktime), 
                 max(plot_data1$soaktime), 
                 length.out = 100),
  waterAct = "west"
)
pred_west$count <- predict(model_west, newdata = pred_west, type = "response")



# Combine predictions
pred_data <- rbind(pred_east, pred_west)


# write the table before plot it 
succOutput[succOutput$phylum==i & succOutput$water=="west",]$slope <- summary(model_west)$coefficients[2,1]
succOutput[succOutput$phylum==i & succOutput$water=="west",]$zslope <- summary(model_west)$coefficients[2,3]
succOutput[succOutput$phylum==i & succOutput$water=="west",]$pslope <- summary(model_west)$coefficients[2,4]



succOutput[succOutput$phylum==i & succOutput$water=="east",]$slope <- summary(model_east)$coefficients[2,1]
succOutput[succOutput$phylum==i & succOutput$water=="east",]$zslope <- summary(model_east)$coefficients[2,3]
succOutput[succOutput$phylum==i & succOutput$water=="east",]$pslope <- summary(model_east)$coefficients[2,4]

tempPlot <- ggplot(plot_data1 %>% filter(phylum==i) %>%
                     filter(treatment=="iBaseline" | treatment=="mBaseline" |
                              treatment=="fBaseline"),
                   aes(y=count, x=soaktime, color=waterAct))+
  geom_point()+
    geom_line(data = pred_data, linewidth = 1, alpha=0.8) +  # Use the corrected pred_data
  theme_classic() +
  labs(title = "Species Richness Succession") +
  scale_color_manual(values = c(   
    "east" = "#56B4E9",   
    "west" = "#CC79A7"
  )) +
  theme_classic()+
  labs(title = i) +
  theme(legend.position = "none",          
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.ticks = element_blank())

assign(paste0("plot.", i), tempPlot) # change the name from tempPlot to paste0 name with assign. 



}


# succOutputLM <- succOutput
# succOutputGLM <- succOutput
succOutputNegBino <- succOutput
succOutputNegBino$pAdj <- p.adjust(succOutputNegBino$pslope, method = "BH")

  write.csv(succOutputNegBino, "decontamByArms/noSED/successionGLM.NB.csv")


# this is for the base model 
plot.porife <- ggplot(plot_data1 %>% filter(phylum=="Porifera") %>%
                        filter(treatment=="iBaseline" | treatment=="mBaseline" |
                                 treatment=="fBaseline"),
                      aes(y=count, x=years, color=waterAct))+
  geom_point()+
  theme_classic()+
  labs(title = "Porifera") +
  scale_color_manual(values = c(   
    "east" = "#56B4E9",   
    "west" = "#CC79A7"))
+
  theme(legend.position = "none",          
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.ticks = element_blank())




plot.allcount
combined_plot2 <- plot_grid(plot.Arthropoda,plot.Annelida,plot.Bacillariophyta,plot.Mollusca,plot.Rhodophyta,plot.Porifera, ncol = 3, nrow=2) # plot the above together 

