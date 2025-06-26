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
library(ggsignif)



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
metaData$years <- "3 years"
metaData[metaData$phase=="initial",]$years <- "1 year"
metaData[metaData$phase=="middle" | metaData$phase=="resistance" ,]$years <- "2 years"

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


## going to plot the resistance community 
# this is the normal group 
# 1 Arthropoda -> pioneer
plot.artho.rt <- ggplot(plot_data1 %>% filter(phylum=="Arthropoda")%>%
                          filter(phase=="resistance" | phase == "middle"),
                        aes(y=count, x=waterAct, color = waterAct))+
  geom_boxplot()+
  geom_jitter(width = 0.15,size=3, alpha=0.8,aes(colour = treatment))+
  scale_color_manual(values = c("#CC79A7","#56B4E9","#0072B2","#E69F00","#D55E00")) +
  theme_classic()+
  labs(title = "Arthropoda") +
  theme(legend.position = "none",          
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.ticks = element_blank(),
        axis.text.x = element_blank())


# 2 Annelida -> pioneer
plot.anne.rt <- ggplot(plot_data1 %>% filter(phylum=="Annelida")%>%
                      filter(phase=="resistance" | phase == "middle"),
                      aes(y=count, x=waterAct, color = waterAct))+
  geom_boxplot()+
  geom_jitter(width = 0.15,size=3, alpha=0.8,aes(colour = treatment))+
  scale_color_manual(values = c("#CC79A7","#56B4E9","#0072B2","#E69F00","#D55E00")) +
  theme_classic()+
  labs(title = "Annelida") +
  theme(legend.position = "none",          
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.ticks = element_blank(),
        axis.text.x = element_blank())

# 3 Bacillariophyta -> indicator # 3 Bacillariophyta ->waterAct indicator 
plot.bacill.rt <- ggplot(plot_data1 %>% filter(phylum=="Bacillariophyta")%>%
                           filter(phase=="resistance" | phase == "middle"),
                         aes(y=count, x=waterAct, color = waterAct))+
  geom_boxplot()+
  geom_jitter(width = 0.15,size=3, alpha=0.8,aes(colour = treatment))+
  scale_color_manual(values = c("#CC79A7","#56B4E9","#0072B2","#E69F00","#D55E00")) +
  theme_classic()+
  labs(title = "Bacillariophyta") +
  theme(legend.position = "none",          
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.ticks = element_blank(),
        axis.text.x = element_blank())

# 4 Rhodophyta -> indicator
plot.rhodop.rt <- ggplot(plot_data1 %>% filter(phylum=="Rhodophyta")%>%
                           filter(phase=="resistance" | phase == "middle"),
                         aes(y=count, x=waterAct, color = waterAct))+
  geom_boxplot()+
  geom_jitter(width = 0.15,size=3, alpha=0.8,aes(colour = treatment))+
  scale_color_manual(values = c("#CC79A7","#56B4E9","#0072B2","#E69F00","#D55E00")) +
  theme_classic()+
  labs(title = "Rhodophyta") +
  theme(legend.position = "none",          
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.ticks = element_blank())


# 5 Mollusca -> pioneer
plot.mollus.rt <- ggplot(plot_data1 %>% filter(phylum=="Mollusca")%>%
                           filter(phase=="resistance" | phase == "middle"),
                         aes(y=count, x=waterAct, color = waterAct))+
  geom_boxplot()+
  geom_jitter(width = 0.15,size=3, alpha=0.8,aes(colour = treatment))+
  scale_color_manual(values = c("#CC79A7","#56B4E9","#0072B2","#E69F00","#D55E00")) +
  theme_classic()+
  labs(title = "Mollusca") +
  theme(legend.position = "none",          
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.ticks = element_blank())




# 6 Porifera -> indicator
plot.porife.rt <- ggplot(plot_data1 %>% filter(phylum=="Porifera")%>%
                           filter(phase=="resistance" | phase == "middle"),
                         aes(y=count, x=waterAct, color = waterAct))+
  geom_boxplot()+
  geom_jitter(width = 0.15,size=3, alpha=0.8,aes(colour = treatment))+
  scale_color_manual(values = c("#CC79A7","#56B4E9","#0072B2","#E69F00","#D55E00")) +
  theme_classic()

+
  labs(title = "Porifera") +
  theme(legend.position = "none",          
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        axis.ticks = element_blank())+
  scale_x_discrete(labels = c("West", "East"))


combined_plot2.rt <- plot_grid(plot.artho.rt, plot.anne.rt, plot.bacill.rt,
                               plot.mollus.rt,plot.rhodop.rt, plot.porife.rt, ncol = 3, nrow=2) # plot the above together 

