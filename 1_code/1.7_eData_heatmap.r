---------------------------------------------------------------------------------------------------------------
# Figure 1d heatmap 
---------------------------------------------------------------------------------------------------------------

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



---------------------------------------------------------------------------------------------------------------
# Table 1 Mean environmental data by site 
---------------------------------------------------------------------------------------------------------------
## load lib
library(fmsb)
library(tidyverse)
library(vegan)
library(ggord)
library(pairwiseAdonis)

## set WD 
setwd("/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#1 biodiveristy/#2 CRF EXP1/#2 dataAnalysis_all/1_data/environmental")
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
site <- c("PC","LM","CDA", "SK","CI","TPC")
phase <- c("initial", "middle","final")

PCrow <- rbind(edata.d2[(edata.d2$Station=="SM10" & edata.d2$Depth =="Bottom Water"),],
               edata.d2[(edata.d2$Station=="SM9" & edata.d2$Depth =="Middle Water"),],
               edata.d2[(edata.d2$Station=="SM11" & edata.d2$Depth =="Middle Water"),])

LMrow <- rbind(edata.d2[(edata.d2$Station=="SM3" & edata.d2$Depth =="Bottom Water"),],
               edata.d2[(edata.d2$Station=="SM4" & edata.d2$Depth =="Middle Water"),])

CDArow <- rbind(edata.d2[(edata.d2$Station=="SM1" & edata.d2$Depth =="Middle Water"),],
               edata.d2[(edata.d2$Station=="MM8" & edata.d2$Depth =="Surface Water"),])

SKrow <- edata.d2[(edata.d2$Station=="PM4" & edata.d2$Depth =="Bottom Water"),]

CIrow <- edata.d2[(edata.d2$Station=="TM4" & edata.d2$Depth =="Bottom Water"),]

TPCrow <- edata.d2[(edata.d2$Station=="MM5" & edata.d2$Depth =="Surface Water"),]

# rework PCrow to mean data of the same period from different site into one 
# so all 6 sites have same same data points 

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
edata.d3 <- rbind(nPCrow,nLMrow,nCDArow,SKrow,CIrow,TPCrow) # 67 rows in total, all 6 sites have 11:12 rows

# write in richness 
metadata <- metadata.d[,1:5]
names(metadata) <- c("ARMS","OTUs","uniOTUs","perUni","site")

edata.d3[,20:21] <- 0
colnames(edata.d3)[20:21] <- c("rich","ARMS")

edata.d4 <- rbind(edata.d3,edata.d3) # 2 ARMS per site

for (i in 1:67) {
  for (j in 1:2) {
    
    edata.d4[i+(j-1)*67,20] <- metadata%>%filter(site==edata.d3[i,1])%>%.[j,2]
    edata.d4[i+(j-1)*67,21] <- metadata%>%filter(site==edata.d3[i,1])%>%.[j,1]
    

  }
  
}

# input treatment and region in edata.d4
edata.d4$treatment <- "MP"
edata.d4$region <- "west"

edata.d4 <- edata.d4 %>% mutate(region= ifelse(site %in% c("TPC","NP","CI"),"east", region))
edata.d4 <- edata.d4 %>% mutate(treatment= ifelse(site %in% c("LM","SK"),"aquaculture", treatment))
edata.d4 <- edata.d4 %>% mutate(treatment= ifelse(site %in% c("PC","CI"),"sewage", treatment))

edata.d4$treatment <- factor(edata.d4$treatment, levels =  c("MP", "aquaculture","sewage"))
edata.d4$site <- factor(edata.d4$site, levels =  c("TPC", "CDA","SK","LM","NP","CC","CI","PC"))

# now let's see if we can build the same model 
model1 <- lm(rich~Total.Inorganic.Nitrogen..mg.L.+region+treatment,data=edata.d4)
summary(model1)

plot1 <- ggplot(edata.d4, aes(x = Total.Inorganic.Nitrogen..mg.L., y = rich, color = treatment, shape = region,group =1)) +
  geom_point(size = 3) +  # Adjust size of points as needed
  labs(x = "TN", y = "OTUs", title = "TN vs OTU") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal()  

model2 <- lm(rich~Chlorophyll.a..μg.L.,data=edata.d4)
summary(model2)

plot2 <- ggplot(edata.d4, aes(x = Chlorophyll.a..μg.L., y = rich, color = treatment, shape = region,group =1)) +
  geom_point(size = 3) +  # Adjust size of points as needed
  labs(x = "col", y = "OTUs", title = "col vs OTU") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_minimal()  


# run a forloop the get all the rich vs eData out 
# rearrange the alpha data to make it easier to loop 
statsDF <- as.data.frame(matrix(0, nrow= 2, ncol= 13))
row.names(statsDF) <- c("pValue","AdjR2")
names(statsDF) <- names(edata.d4)[7:19]
edata.d5 <- edata.d4[1:67,]


# 1,2,4 it's okay no need to transfer 
# 3 do log(max salinity + 1 - salinity)
# 5: 15 do log(x +1)

for (i in c(1,2)) {
  modelTemp <- lm(rich~ . ,edata.d4[,c("rich",names(statsDF)[i])])
  statsDF[1,i] <- summary(modelTemp)$coefficients[2,4] 
  statsDF[2,i] <- summary(modelTemp)$adj.r.squared
  
}


for (i in c(3:13)) {
  formula_str <- as.formula(paste("log(",names(statsDF)[i], "+1)~ rich"))
  modelTemp <- lm(formula_str, edata.d5)
  statsDF[1,i] <- summary(modelTemp)$coefficients[2,4] 
  statsDF[2,i] <- summary(modelTemp)$adj.r.squared
}  

statsLong <- as.data.frame(t(statsDF))

statsLongSorted <- statsLong %>%
  arrange(.[[1]]) 



# from the above result we can see after puting all the data point out the model is not able to predict 
# as good as it were when just using the mean 


# are treatment differetn from each other in some parameter? 

# check normality 
hist(edata.d5[,7]) # 1 DO normal
hist(edata.d5[,8]) # 2 DO saturation normal
hist(edata.d5[,9]) # 5 tubidity 
hist(edata.d5[,10]) # 6 suspended solids right skew 
hist(edata.d5[,11]) # 7 TIN right skew right skew
hist(edata.d5[,12]) # 8 Ammonia right skew
hist(edata.d5[,13]) # 9 Nitrite right skew
hist(edata.d5[,14]) # 10 Nitrate right skew
hist(edata.d5[,15]) # 11 P right skew
hist(edata.d5[,16]) # 12 E.coli right skew
hist(edata.d5[,17]) # 13 Faecal Coli right skew
hist(edata.d5[,18]) # 14 Chlorophyll a right skew
hist(edata.d5[,19]) # 15 Phaeo right skew


# check normality transformed 
hist(edata.d5[,7]) # 1 DO normal
hist(edata.d5[,8]) # 2 DO saturation normal
hist(log(edata.d5[,9])) # 5 tubidity 
hist(log(edata.d5[,10])) # 6 suspended solids, normal after transfer 
hist(log(edata.d5[,11])) # 7 TIN right skew right skew before log transfer but also weird afterwards 
hist(log(edata.d5[,12]))  # 8 Ammonia right skew
hist(log(edata.d5[,13]))  # 9 Nitrite right skew
hist(log(edata.d5[,14]))  # 10 Nitrate right skew
hist(log(edata.d5[,15]))  # 11 P right skew
hist(log(edata.d5[,16])) # 12 E.coli right skew
hist(log(edata.d5[,17])) # 13 Faecal Coli right skew
hist(log(edata.d5[,18])) # 14 Chlorophyll a right skew
hist(log(edata.d5[,19])) # 15 Phaeo right skew

# not all perfect but largly looks better 
# 1,2, it's okay no need to transfer 
# 5: 15 do log(x +1)

edata.d5$n2p <- edata.d5$Total.Inorganic.Nitrogen..mg.L./edata.d5$Orthophosphate.Phosphorus..mg.L.
edata.d5[edata.d5$n2p<16,]$site

plot3 <- ggplot(edata.d5, aes(x = site, y = Total.Inorganic.Nitrogen..mg.L., color = treatment, shape = region)) +
  geom_boxplot() +  # Adjust size of points as needed
  labs(x = "site", y = "TIN") +
  theme_minimal()  

plot3.a <- ggplot(edata.d5, aes(x = site, y = n2p, color = treatment, shape = region)) +
  geom_boxplot() +  # Adjust size of points as needed
  labs(x = "site", y = "TIN") +
  theme_minimal()  


model3.1 <-  lm(Orthophosphate.Phosphorus..mg.L.~treatment+region,data=edata.d5)
summary(model3.1)


plot3.1 <- ggplot(edata.d5, aes(x = site, y = Orthophosphate.Phosphorus..mg.L., color = treatment, shape = region)) +
  geom_boxplot() +  # Adjust size of points as needed
  labs(x = "TP", y = "site") +
  theme_minimal()  



model4 <-  lm(Chlorophyll.a..μg.L.~treatment+region,data=edata.d5)
summary(model4)

plot4 <- ggplot(edata.d5, aes(x = site, y = Chlorophyll.a..μg.L., color = treatment, shape = region)) +
  geom_boxplot() +  # Adjust size of points as needed
  labs(x = "col", y = "sites") +
  theme_minimal()  


model4.1 <-  lm(Phaeo.pigments..μg.L.~treatment+region,data=edata.d5)
summary(model4.1)

plot4.1 <- ggplot(edata.d5, aes(x = site, y = Phaeo.pigments..μg.L., color = treatment, shape = region)) +
  geom_boxplot() +  # Adjust size of points as needed
  labs(x = "phaeo", y = "sites") +
  theme_minimal()  


model5 <-  lm(E..coli..cfu.100mL.~treatment+region,data=edata.d5)
summary(model5)

plot5 <- ggplot(edata.d5, aes(x = site, y = E..coli..cfu.100mL., color = treatment, shape = region)) +
  geom_boxplot() +  # Adjust size of points as needed
  labs(x = "eco", y = "sites") +
  theme_minimal()  

model5.1 <-  lm(Faecal.Coliforms..cfu.100mL.~treatment+region,data=edata.d5)
summary(model5.1)

plot5.1 <- ggplot(edata.d5, aes(x = site, y = Faecal.Coliforms..cfu.100mL., color = treatment, shape = region)) +
  geom_boxplot() +  # Adjust size of points as needed
  labs(x = "Faecal", y = "sites") +
  theme_minimal()  

# can i build a for loop again to plot out all the stats? 


statsTreReg <- as.data.frame(matrix(0, nrow= 3, ncol= 13))
row.names(statsTreReg) <- c("pA","pSew","pRegion")
names(statsTreReg) <- names(edata.d4)[7:19]

# 1,2,4 it's okay no need to transfer 
# 3 do log(max salinity + 1 - salinity)
# 5: 15 do log(x +1)

for (i in c(1,2)) {
  formula_str <- as.formula(paste(names(statsTreReg)[i], "~ treatment + region"))
  modelTemp <- lm(formula_str, edata.d5)
  statsTreReg[1,i] <- summary(modelTemp)$coefficients[2,4] 
  statsTreReg[2,i] <- summary(modelTemp)$coefficients[3,4] 
  statsTreReg[3,i] <- summary(modelTemp)$coefficients[4,4] 

}

for (i in c(3:13)) {
  formula_str <- as.formula(paste("log(",names(statsTreReg)[i], "+1)~ treatment + region"))
  modelTemp <- lm(formula_str, edata.d5)
  statsTreReg[1,i] <- summary(modelTemp)$coefficients[2,4] 
  statsTreReg[2,i] <- summary(modelTemp)$coefficients[3,4] 
  statsTreReg[3,i] <- summary(modelTemp)$coefficients[4,4] 
  
}


statsTreRegLong <- as.data.frame(t(statsTreReg))
statsTreRegLong$adjPvalueBH.A <- p.adjust(statsTreRegLong$pA, method = "BH")
statsTreRegLong$adjPvalueBH.Sew <- p.adjust(statsTreRegLong$pSew, method = "BH")
statsTreRegLong$adjPvalueBH.Region <- p.adjust(statsTreRegLong$pRegion, method = "BH")



statsTreRegLongSort <- statsTreRegLong %>%
  arrange(.[[6]]) 


# Define c1 as a vector of colors
c1 <- c("#1F77B4", "#17BECF","#D62728","#FF7F0E",  "#E377C2",  "#9467BD", "#8C564B","#2CA02C")


set.seed(123)
ord<-metaMDS(edata.d5[,c(7:19)])
plotx1 <- ggord(ord,edata.d5$site, size =2, alpha_el = 0.1, poly = FALSE)+
  scale_color_manual(values=c(c1[1:8]),limits=c("TPC","CDA","SK","LM","NP","CC","CI","PC"))+
  scale_fill_manual(values=c(c1[1:8]),limits=c("TPC","CDA","SK","LM","NP","CC","CI","PC"))+
  coord_cartesian(xlim = c(-0.8,0.5), ylim = c(-0.4,0.4))+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

plotx2 <- ggord(ord,edata.d5$region, size =2, alpha_el = 0.1, poly = FALSE)+
  scale_color_manual(values=c(c1[1],c1[3]),limits=c("west", "east"))+
  scale_fill_manual(values=c(c1[1],c1[3]),limits=c("west", "east"))+
  coord_cartesian(xlim = c(-0.8,0.5), ylim = c(-0.4,0.4))+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

plotx3 <- ggord(ord,edata.d5$treatment, size =2, alpha_el = 0.1, poly = FALSE)+
  scale_color_manual(values=c(c1[1:4]),limits=c("sewage", "sedimentation","aquaculture","MP"))+
  scale_fill_manual(values=c(c1[1:4]),limits=c("sewage", "sedimentation","aquaculture","MP"))+
  coord_cartesian(xlim = c(-0.8,0.5), ylim = c(-0.4,0.4))+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


env.mat<-edata.d5[,c(7:19)]
env.dist<-vegdist(env.mat, method="euclidian")



set.seed(100)
modelx2 <- adonis2(env.dist~treatment+region,data=edata.d5, permutations = 999, strata = edata.d5$Dates,  by = "term" )
# Df SumOfSqs      R2      F Pr(>F)   
# treatment         3  8531100 0.13704 4.6603  0.003 **
# treatment:region  4  3076319 0.04942 1.2604  0.265   
# Residual         83 50646194 0.81355                 
# Total            90 62253614 1.00000                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


set.seed(100)
modelx2.1 <- pairwise.adonis(x=edata.d5[,c(7:19)],factors=edata.d5$treatment,
                              sim.function='vegdist', sim.method='euclidian',p.adjust.m='bonferroni')

#result: 
# pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 sewage vs aquaculture  1  739261.4 0.6068625 0.01360469   0.460      1.000    
# 2          sewage vs MP  1 5802753.2 7.1567943 0.14559115   0.001      0.003   *
# 3     aquaculture vs MP  1 2467000.8 5.2979830 0.11201287   0.009      0.027   .





### ok let's model turbidity and TSS and see what's going on
TUBvsTSS <- ggplot(edata.d5, aes(x=Turbidity..NTU., y=Suspended.Solids..mg.L., group=treatment)) +
  scale_shape_manual(values = c(0,1,2,3)) +
  geom_point(aes(shape=treatment))+
  stat_ellipse(aes(fill = treatment), type = "norm", level = 0.99, geom = "polygon", alpha = 0.2) +  # Add ellipses
  theme_minimal() 


TUBvsTSS_OUTremoved <- ggplot(edata.d5%>%filter(Turbidity..NTU.<50), aes(x=Turbidity..NTU., y=Suspended.Solids..mg.L., group=treatment)) +
  scale_shape_manual(values = c(0,1,2,3)) +
  geom_point(aes(shape=treatment))+ø
  stat_ellipse(aes(fill = treatment), type = "norm", level = 0.99, geom = "polygon", alpha = 0.2) +  # Add ellipses
  theme_minimal() 

# see if any group stands out in the adonis model 
set.seed(100)
adonis_result <- adonis2(edata.d5[,c(11:12)] ~ treatment, data=edata.d5, method = "bray")
adonis_result.pw <- pairwise.adonis2(edata.d5[,c(11:12)] ~ treatment, data=edata.d5, by = "terms")


## by the end write everything out 
write.csv(statsTreRegLongSort, "NoSed/edataDifference.csv",row.names = T)
write.csv(statsLongSorted, "NoSed/AlphaEdataBreak.csv",row.names = T)


