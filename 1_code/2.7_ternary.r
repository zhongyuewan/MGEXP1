## set WD 
#mac 
setwd("/Users/moicomputer/Library/CloudStorage/OneDrive-TheUniversityofHongKong-Connect/#2 PhD/#2 research project/#1 biodiveristy/#2 CRF EXP1/#2 dataAnalysis_all/1_data")
rm(list=ls())

## load lib
library(dplyr)
library(ggplot2)
library(vegan)
library(ggtern)
library(cowplot)

## get data
freqTable <- read.csv('decontamByArms/noSED/freqTable.csv', header = TRUE, row.names = 1)
metadata <- read.csv('decontamByArms/noSED/metadata.csv', row.names = 1)

## fix table names and stuffs 
names(freqTable) <- gsub("X","", names(freqTable))

# build distance matrix 
Bcc_matrix <-as.matrix(vegdist(t(freqTable), method = "jaccard")) # bray with sr transfered data 

## jaccard pairwise is Bcc, now i need to calculate Brich
Brich_matrix <- matrix(0, nrow = ncol(freqTable), ncol = ncol(freqTable))
colnames(Brich_matrix) <- c(names(freqTable))
row.names(Brich_matrix)<- names(freqTable)

for (i in 1:(ncol(freqTable) - 1)) {
  for (j in (i + 1):ncol(freqTable)) {
    all_reads <- sum(freqTable[, i] | freqTable[, j])  # Number of shared reads
    unique_reads_A <- sum(freqTable[, i] & !freqTable[, j])  # Number of unique reads in sample A
    unique_reads_B <- sum(!freqTable[, i] & freqTable[, j])  # Number of unique reads in sample B
    
    # Calculate the pairwise value using the formula you provided
    pairwise_value <- abs((unique_reads_A - unique_reads_B) / all_reads)
    
    # Store the pairwise value in the matrix
    Brich_matrix[i, j] <- pairwise_value
    Brich_matrix[j, i] <- pairwise_value
  }
}

# makes sure names are the same before 
all(names(Bcc_matrix) == names(Brich_matrix))
all(row.names(Bcc_matrix) == row.names(Brich_matrix))


Brp_matrix <- Bcc_matrix - Brich_matrix


## create different subset for 
# west vs east
westArms <- row.names(metadata[metadata$site1=="CDA",])
eastArms <- row.names(metadata[metadata$site1=="TPC",])

# initial / middle / final 
westArms.ib <- row.names(metadata[metadata$site1=="CDA" & metadata$phase=="initial", ])
eastArms.ib <- row.names(metadata[metadata$site1=="TPC" & metadata$phase=="initial", ])
westArms.mb <- row.names(metadata[metadata$site1=="CDA" & metadata$phase=="middle", ])
eastArms.mb <- row.names(metadata[metadata$site1=="TPC" & metadata$phase=="middle", ])
westArms.fb <- row.names(metadata[metadata$site1=="CDA" & metadata$phase=="final", ])
eastArms.fb <- row.names(metadata[metadata$site1=="TPC" & metadata$phase=="final", ])



# resistance 
westArms.RT <- row.names(metadata[metadata$site1=="CDA" & metadata$phase=="resistance",])
eastArms.RT <- row.names(metadata[metadata$site1=="TPC" & metadata$phase=="resistance",])


# resilience 
westArms.RL <- row.names(metadata[metadata$site1=="CDA" & metadata$phase=="resilience",])
eastArms.RL <- row.names(metadata[metadata$site1=="TPC" & metadata$phase=="resilience",])

# let's do 1. west middle vs (initial + resistant) 
west.RT <- data.frame(matrix(0, ncol=4, nrow=8))
east.RT <- data.frame(matrix(0, ncol=4, nrow=8))

colnames(west.RT) <- c("ARMS","Sim", "Brp", "Brich")
west.RT$ARMS <- c(westArms.RT, westArms.RT)
west.RT$site <- metadata[west.RT$ARMS,]$site2
west.RT$pair [1:4] <- "79" 
west.RT$pair [5:8] <- "81" 

colnames(east.RT) <- c("ARMS","Sim", "Brp", "Brich")
east.RT$ARMS <- c(eastArms.RT, eastArms.RT)
east.RT$site <- metadata[east.RT$ARMS,]$site2
east.RT$pair [1:4] <- "69" 
east.RT$pair [5:8] <- "70" 


west.RT$treatment <- metadata[west.RT$ARMS,]$treatment
west.RT$water <- "west"

east.RT$treatment <- metadata[east.RT$ARMS,]$treatment
east.RT$water <- "east"


west.RT[west.RT$pair==79,]$Sim <- 1-Bcc_matrix[west.RT[1:4,]$ARMS,"79"]
west.RT[west.RT$pair==81,]$Sim <- 1-Bcc_matrix[west.RT[5:8,]$ARMS,"81"]

west.RT[west.RT$pair==79,]$Brp <- Brp_matrix[west.RT[1:4,]$ARMS,"79"]
west.RT[west.RT$pair==81,]$Brp <- Brp_matrix[west.RT[5:8,]$ARMS,"81"]

west.RT[west.RT$pair==79,]$Brich <- Brich_matrix[west.RT[1:4,]$ARMS,"79"]
west.RT[west.RT$pair==81,]$Brich <- Brich_matrix[west.RT[5:8,]$ARMS,"81"]

east.RT[east.RT$pair==69,]$Sim <- 1-Bcc_matrix[east.RT[1:4,]$ARMS,"69"]
east.RT[east.RT$pair==70,]$Sim <- 1-Bcc_matrix[east.RT[5:8,]$ARMS,"70"]

east.RT[east.RT$pair==69,]$Brp <- Brp_matrix[east.RT[1:4,]$ARMS,"69"]
east.RT[east.RT$pair==70,]$Brp <- Brp_matrix[east.RT[5:8,]$ARMS,"70"]

east.RT[east.RT$pair==69,]$Brich <- Brich_matrix[east.RT[1:4,]$ARMS,"69"]
east.RT[east.RT$pair==70,]$Brich <- Brich_matrix[east.RT[5:8,]$ARMS,"70"]

RT <- rbind(west.RT, east.RT)


# let's plot RL but this time just do one dataset 
# in resilience, we look at final vs inital/middle/treatments 
RL <- data.frame(matrix(0, ncol=5, nrow=16))
colnames(RL) <- c("ARMS","Sim", "Brp", "Brich", "pair")
RL$ARMS <- c(westArms.RL, westArms.RL,
             eastArms.RL, eastArms.RL)
RL$site <- metadata[RL$ARMS,]$site2
RL$pair[1:4] <- "77" 
RL$pair[5:8] <- "78" 
RL$pair[9:12] <- "71" 
RL$pair[13:16] <- "72" 

RL$water <- "west"
RL$water[9:16] <- "east" 

RL[RL$pair==77,]$Sim <- 1-Bcc_matrix[RL[1:4,]$ARMS,"77"]
RL[RL$pair==78,]$Sim <- 1-Bcc_matrix[RL[5:8,]$ARMS,"78"]
RL[RL$pair==71,]$Sim <- 1-Bcc_matrix[RL[9:12,]$ARMS,"71"]
RL[RL$pair==72,]$Sim <- 1-Bcc_matrix[RL[13:16,]$ARMS,"72"]

RL[RL$pair==77,]$Brp <- Brp_matrix[RL[1:4,]$ARMS,"77"]
RL[RL$pair==78,]$Brp <- Brp_matrix[RL[5:8,]$ARMS,"78"]
RL[RL$pair==71,]$Brp <- Brp_matrix[RL[9:12,]$ARMS,"71"]
RL[RL$pair==72,]$Brp <- Brp_matrix[RL[13:16,]$ARMS,"72"]

RL[RL$pair==77,]$Brich <- Brich_matrix[RL[1:4,]$ARMS,"77"]
RL[RL$pair==78,]$Brich <- Brich_matrix[RL[5:8,]$ARMS,"78"]
RL[RL$pair==71,]$Brich <- Brich_matrix[RL[9:12,]$ARMS,"71"]
RL[RL$pair==72,]$Brich <- Brich_matrix[RL[13:16,]$ARMS,"72"]

RL$treatment <- metadata[RL$ARMS,]$treatment

# plotting
# resistant 
RT[RT$treatment=="eutrophication",]$treatment <- "Sewage"
RT[RT$treatment=="fishFarm",]$treatment <- "Mariculture"

RT$treatment <- factor(RT$treatment, level = c("Sewage", "Mariculture"))
RT.plot <- ggtern(data = RT, aes(Sim, Brp, Brich, color = treatment)) +
  geom_point(
    alpha = 1,
    size = 3,
  ) +
  theme_bw(base_size = 20)+
  labs(
    T = "Replacement",    
    L = "Similarity",     
    R = "Richness" 
  ) +
  scale_color_manual(values = c(   
    "Mariculture" = "#E69F00",   
    "Sewage" = "#D55E00"))+
  theme(
    tern.axis.title.T = element_blank(),
    tern.axis.title.L = element_blank(),
    tern.axis.title.R = element_blank(),
    legend.position="none")


#mid baseline vs resilience  
RL[RL$treatment=="eutrophication",]$treatment <- "Sewage"
RL[RL$treatment=="fishFarm",]$treatment <- "Mariculture"
RL$treatment <- factor(RL$treatment, level = c("Sewage", "Mariculture"))

RL.plot <- ggtern(data = RL, aes(Sim, Brp, Brich, color = treatment)) +
  geom_point(
    alpha = 1,
    size = 3,

  ) +
  theme_bw(base_size = 20)+
  labs(
    T = "Replacement",    
    L = "Similarity",     
    R = "Richness" 
  ) +
  scale_color_manual(values = c(   
    "Mariculture" = "#E69F00",   
    "Sewage" = "#D55E00"))+
  theme(
    tern.axis.title.T = element_blank(),
    tern.axis.title.L = element_blank(),
    tern.axis.title.R = element_blank(),
    legend.position="none")


combined_plot <- plot_grid(RT.plot+  theme(legend.position="none"),RL.plot+  theme(legend.position="none")) # plot the above together 

### run some stats 
RT$phase <- "RT"
RL$phase <- "RL"

RRR <- rbind(RT,RL)

hist(RRR$Brp)
hist(RRR$Sim)
hist(RRR$Brich)

shapiro.test(RRR$Brp)
shapiro.test(RRR$Sim)
shapiro.test(RRR$Brich)

model1 <- t.test(Brp~phase,RRR, var.equal = T)
model1 # RL is significantly more similar than RT

model2 <- t.test(Sim~phase,RRR, var.equal = T)
model2 # RL is significantly more similar than RT

model3 <- t.test(Brich~phase,RRR, var.equal = T)
model3 # RL is significantly more similar than RT

### make visulization 
RRS <- RRR[,-3:-4]
RRB <- RRR[,c(-2,-4)]
RRBr <- RRR[,-2:-3]

RRS$t <- "Sim"
RRB$t <- "Brp"
RRBr$t <- "Brich"

names(RRS)[2] <- "number"
names(RRB)[2] <- "number"
names(RRBr)[2] <- "number"

RRR.plot <- rbind(RRS,RRB,RRBr)


RRR.plot$t <- factor(RRR.plot$t, level=c("Sim","Brp","Brich"))
RRR.plot$phase <- factor(RRR.plot$phase, level = c("RT","RL")) 
plotANOVA <- ggplot(RRR.plot, aes(x = t, y = number, color = phase)) +
  geom_boxplot() +  # Adjust size of points as needed
  labs(x = "Similarity, Replacment, Richness", y = "percentage", title = "Differences in compositional changes") +
  theme_minimal()   






#####################################
# Your data
your_data <- c(0.1310137, 0.1173803, 0.3212716, 0.1037092, 0.1999344, 0.2066021, 0.1248486, 0.1087343, 0.3748278, 0.1348294, 0.2448283, 0.2246577, 0.1181517, 0.1121104,
               0.2886228, 0.1413889, 0.2389521, 0.2108983, 0.1202987, 0.1036302, 0.2626086, 0.1171945, 0.2626390, 0.2454878, 0.2032722, 0.2719520, 0.2281768, 0.2878339,
               0.2621613, 0.2776882, 0.2932198, 0.2934927, 0.2507561, 0.2407569, 0.3909133, 0.3293421, 0.1892459, 0.2398683, 0.2728218, 0.2571020, 0.2673307, 0.2573093,
               0.2783222, 0.2942617, 0.3082886, 0.2445311, 0.3429616, 0.2585495)

# Example independent variable (replace with your actual data)
x <- 1:length(your_data)

# Check distribution
hist(your_data, breaks = 20, main = "Histogram of Data", xlab = "Values")
qqnorm(your_data)
qqline(your_data)

# Apply a transformation if necessary (e.g., log transformation)
transformed_data <- log(your_data + 1)

# Fit a linear model
model <- lm(your_data ~ x)

# View the model summary
summary(model)

# Check assumptions
shapiro.test(residuals(model))  # Normality of residuals
plot(model, which = 1)          # Residuals vs Fitted plot
plot(x, transformed_data)       # Scatterplot of x vs y
abline(model, col = "red")      # Add regression line



model2 <- lm(Brp~phase,RRR)
summary(model2) # RT has significantly more replacement than RL 

model3 <- lm(Brich~phase,RRR)
summary(model3) # no significant difference on richness, so community difference were driven by replacement  


