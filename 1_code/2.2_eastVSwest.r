## load lib
library(tidyverse)
library(factoextra)
library(RColorBrewer)

## set WD 
# pc 
setwd("E:/OneDrive - connect.hku.hk/#2 PhD/#2 research project/#1 biodiveristy/#2 CRF EXP1/#2 dataAnalysis_all/1_data/environmental/")
# mac
setwd("/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#1 biodiveristy/#2 CRF EXP1/#2 dataAnalysis_all/1_data/environmental")
rm(list=ls())

# read file 
edata <- read.csv("edata.csv", header = T)

names(edata) # check what i need 
# General -> Dissolved.Oxygen..mg[26], Dissolved.Oxygen...saturation[27]
# General -> Salinity[18], pH[20]
# sedimentation ->  Turbidity [9], Suspended.Solids[15]
# aquaculture -> Total.Inorganic.Nitrogen[13],Ammonia.Nitrogen[29],
# aquaculture -> Nitrite.Nitrogen[22],Nitrate.Nitrogen[23],Orthophosphate.Phosphorus[19],
# sewage -> E..coli..cfu[25], Faecal.Coliforms[24]
# sewage -> Chlorophyll.a..[28], Phaeo.pigments [19]

# also check Silica for bacillariophyta 

edata.d1 <- edata[,c(31,32,3,4,5,6,26,27,9,15,13,29,22,23,19,25,24,28,21)]

siteP <- read.csv("sitePofile.csv")

# filter dates 
edata.d2 <- edata.d1 
  
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


######################################### below are old codes #############
####from those code, i did not mean the temporal/spacial replicate first ##  
#### which i should, otherwise outliner will be more outliner #############

############## but then this 

############## check outliner 
edata.d5 %>% filter(site=="TPC") %>% .$Turbidity..NTU. %>% mean() # 13.77 mean everything
edata.d5 %>% filter(site=="TPC") %>%
  group_by(Dates) %>%
  summarise(mean_value = mean(Turbidity..NTU.)) %>% .$mean_value %>% mean
##############################
# so acutally it does not matter


# build a new df to store the mean and sd 
edata.d6 <-   data.frame(matrix(0, nrow=10, ncol=13))
row.names(edata.d6) <-  c(unique(edata.d5$site), "mean", "sd")
names(edata.d6) <- names(edata.d5)[7:19]

# make a function to mean all the things 
meanEdata <- function (df){
  for (i in 1:8){
for (j in 1:13){
  edata.d6[i,j] <- mean(df[df$site == row.names(edata.d6[i,]),][,j+6])
}
  }
  return(edata.d6)
}


edata.d6 <- meanEdata(edata.d5)

# calculate group mean and sd
for (i in 1:13){ #mean
edata.d6[9,i] <- mean(edata.d6[1:8,i])

}

for (i in 1:13){ #sd
  edata.d6[10,i] <- sd(edata.d6[1:8,i])
  
}


# turn this table into -mena / sd 
edata.d7 <- edata.d6
for (i in 1:8){ #how many of sd
  for (j in 1:13){
  edata.d7[i,j] <- (edata.d6[i,j] - edata.d6[9,j])/edata.d6[10,j]
}
}

# add treatment/water  
edata.d7$water <- "west"
edata.d7$water[row.names(edata.d7) == "NP" | row.names(edata.d7) == "TPC" |
                 row.names(edata.d7) == "CI" |row.names(edata.d7) == "SK"] <- "east"

edata.d7$treatment <- "baseline"
edata.d7$treatment[row.names(edata.d7) == "NP"|row.names(edata.d7) == "CC"] <- "sedimentation"
edata.d7$treatment[row.names(edata.d7) == "PC"|row.names(edata.d7) == "CI"] <- "eutrophication"
edata.d7$treatment[row.names(edata.d7) == "LM"|row.names(edata.d7) == "SK"] <- "fishfarm"

edata.d7$water[9:10] <- "meansd"
edata.d7$treatment[9:10] <- "meansd"

edata.d7 <- mutate(edata.d7, treatmentWater = paste0(water,treatment))

# make clear df to plot figure
edata.d8 <- edata.d7[1:8,] # remove mean/sd
edata.fi <- edata.d8
edata.fi[,c(1,2,3,4)] <- -edata.fi[,c(1,2,3,4)] # reverse 1,2,3 which is o2, salinity, pH 

# reorder them by east/west & impact 
customOrder <- c("CDA","CC","LM","PC","TPC","NP","SK","CI")
edata.plot <- edata.fi[match(customOrder,rownames(edata.fi)),]
edata.plot$Mix_impact <- rowMeans(edata.plot[13])
write.csv(edata.plot,"edata.plot_3year.csv", row.names = T)
write.csv(edata.d7,"edata.plotBeforeNormal_3year.csv", row.names = T)

edata.d8 <- edata.d7[c(2,3),]


# heat map 
# show color 
myColor <- colorRampPalette(c("deepskyblue", "white", "coral3"))(80)


# plot
plot1 <- heatmap(as.matrix(t(edata.plot[,c(1:13,17)])), Colv = NA, Rowv = NA, 
                 scale="row", col = myColor) 

# it's not too much different from the old one but now i can see better trends 

## what if i want to do some light stat just to see if any site is really 

# reorder them by impact & east/west  
customOrder <- c("TPC","CDA", "NP","CC","SK","LM","CI","PC")
edata.plot <- edata.fi[match(customOrder,rownames(edata.fi)),]
edata.plot$Mix_impact <- rowMeans(edata.plot[,1:13])

plot2 <- heatmap(as.matrix(t(edata.plot[,c(1:13,17)])), Colv = NA, Rowv = NA, 
                 scale="row", col = myColor) 


# check stats and see if TPC is better than CDA in ChlA and nitorgen 
model1 <- lm(Nitrite.Nitrogen..mg.L. ~site+phase, data = edata.d5 
             %>% filter(site == "NP" | site == "TPC"))
summary(model1)

plot3 <- ggplot(data = edata.d5 %>% filter(site == "CDA" | site == "TPC"), 
                aes(y = Total.Inorganic.Nitrogen..mg.L. , x = site, color = phase)) +
  geom_boxplot() +
  labs(x = "site", y = "ChlA", color = "Depth") 


### check total inortganic nitrogen temporal changes over 3 years in CDA and TPC 

# filter dates 


# filter depth, 
edata3y.d1 <- edata.d1 %>%
  mutate(mergecell=paste0(Station,Depth)) %>% 
  filter(mergecell %in% siteP$mergecell)

edata3y.d2 <- edata3y.d1
for (i in 7:19){# remove "<"
  edata3y.d2[,i] <- as.numeric(ifelse(grepl("<", edata3y.d1[,i]), 0, edata3y.d1[,i]))
}

# mutate a unique id for each sampling point 
edata3y.d3 <- mutate(edata3y.d2, ID = paste0(site, row_number()))
row.names(edata3y.d3) <- edata3y.d3$ID



# only CDA/TPC
edata3y.d4 <- edata3y.d3 %>% filter (site == "CDA" | site == "TPC")

model3y.1 <- lm(Total.Inorganic.Nitrogen..mg.L.~site, data=edata3y.d4)
summary(model3y.1) # p = 0.000956

model3y.2 <- lm(Chlorophyll.a..μg.L.~site+phase, data=edata3y.d4)
summary(model3y.2) # p = 0.0466

model3y.3 <- lm(Orthophosphate.Phosphorus..mg.L.~site+phase, data=edata3y.d4)
summary(model3y.3) # p = 0.00618

model3y.4 <- lm(E..coli..cfu.100mL.~site+phase, data=edata3y.d4)
summary(model3y.4) # p = 0.439 not significant 



plot3y.1 <- ggplot(data = edata3y.d4, 
                aes(y = Total.Inorganic.Nitrogen..mg.L. , x = site, color = phase)) +
  geom_boxplot() +
  labs(x = "site", y = "TIN", color = "site") +
  theme_classic()
  
  
plot3y.2 <- ggplot(data = edata3y.d4, 
                     aes(y = Chlorophyll.a..μg.L. , x = site, color = phase )) +
  geom_boxplot() +
  labs(x = "site", y = "ChlA", color = "site") +
  theme_classic()


## what's significant in 3 years?
# run a forloop the get all the rich vs eData out 
# rearrange the alpha data to make it easier to loop 
statsDF <- as.data.frame(matrix(0, nrow= 4, ncol= 13))
row.names(statsDF) <- c("CDA","TPC","pValue","AdjR2")
names(statsDF) <- names(edata.plot)[1:13]
edata.d5 <- edata3y.d4


# 1,2 it's okay no need to transfer 
# 3: 13 do log(x +1)

for (i in c(1,2)) {
  formula_str <- as.formula(paste(names(statsDF)[i], "~ site"))
  modelTemp <- lm(formula_str, edata.d5)
  statsDF[3,i] <- summary(modelTemp)$coefficients[2,4] 
  statsDF[4,i] <- summary(modelTemp)$adj.r.squared
  statsDF[1,i] <- edata.d5%>%filter(site=="CDA") %>% .[,6+i] %>% mean(.) 
  statsDF[2,i] <- edata.d5%>%filter(site=="TPC") %>% .[,6+i] %>% mean(.)
  
}

for (i in c(3:13)) {
  formula_str <- as.formula(paste("log(",names(statsDF)[i], "+1)~ site"))
  modelTemp <- lm(formula_str, edata.d5)
  statsDF[3,i] <- summary(modelTemp)$coefficients[2,4] 
  statsDF[4,i] <- summary(modelTemp)$adj.r.squared
  statsDF[1,i] <- edata.d5%>%filter(site=="CDA") %>% .[,6+i] %>% mean(.) 
  statsDF[2,i] <- edata.d5%>%filter(site=="TPC") %>% .[,6+i] %>% mean(.)
  
  
  }  

statsLong <- as.data.frame(t(statsDF))
statsLong$adjPvalueBH <- p.adjust(statsLong$pValue, method = "BH")



statsLongSorted <- statsLong %>%
  arrange(.[[5]]) 



write.csv(statsLongSorted, "../decontamByArms/noSED/AlphaEdataBreak3year.csv",row.names = T)


