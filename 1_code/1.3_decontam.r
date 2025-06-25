
## load lib
library(Biostrings)
library(dplyr)
library(biomformat)


# mac
setwd("/filepath")
rm(list=ls())

## get data
# convert biom to csv 
OTU_reads <- read_biom(r"(feature-table.biom)")
OTU_table <- as.data.frame(as.matrix(biom_data(OTU_reads)))
freqTable <- OTU_table

metaData <- read.table('sample-metadata.tsv', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, comment.char = '')
sequenceTable <- readDNAStringSet('dna-sequences.fasta')


## play with data
# 1. find columns that contain NTC in header
# 2. sum rows that contain NEG in header to column NEGSUM, 
# 3. drop individual NEG columns
## negColumns <- grep('COTR', names(freqTable), value = T)
## freqTable$NEGSUM <- rowSums(freqTable[negColumns])
## freqTable <- freqTable[, !names(freqTable) %in% negColumns]

# there are four control samples that controls different samples 
# first is to break the big dataset into four different subset 
nc1Rows <- row.names(metaData)[metaData[,4]=="nc1"]
nc2Rows <- row.names(metaData)[metaData[,4]=="nc2"]
nc3Rows <- row.names(metaData)[metaData[,4]=="nc3"]
nc4Rows <- row.names(metaData)[metaData[,4]=="nc4"]

# Look at how many features are picked up by the ngative control
# 1. create list of index names (ZOTU numbers) for which freqTable$NEGSUM > 0, 
# 2. print taxonomy of ZOTUs in list
negZOTUs1 <- rownames(freqTable)[which(freqTable[nc1Rows]$nc1 > 0)]
length(negZOTUs1)
negZOTUs2 <- rownames(freqTable)[which(freqTable[nc2Rows]$nc2 > 0)]
length(negZOTUs2)
negZOTUs3 <- rownames(freqTable)[which(freqTable[nc3Rows]$nc3 > 0)]
length(negZOTUs3)
negZOTUs4 <- rownames(freqTable)[which(freqTable[nc4Rows]$nc4 > 0)]
length(negZOTUs4)


# identify how many negative features are picked by in the positive samples. 
#1 nc1
negTotal1 <- data.frame(matrix(nrow=length(negZOTUs1),ncol=5))
rownames(negTotal1) <- negZOTUs1
colnames(negTotal1) <- c("readsInNeg","%ofTotalRead", "meanReads", "PresentInPositive", "%ofPositivepresent")

for (x in negZOTUs1) {
  if (x %in% negZOTUs1) {
    negValue <- freqTable[x, 'nc1']
    negValuePerc <- negValue / rowSums(freqTable[x,nc1Rows]) * 100
    rowMeanValue <- rowMeans(freqTable[x,nc1Rows])
    positiveDetections <- sum(freqTable[x, nc1Rows] > 0)-1
    positiveDetectionsPerc <- (sum(freqTable[x,nc1Rows] > 0)-1) / (length(nc1Rows)-1) * 100
    negTotal1[x,] <- c(negValue, negValuePerc, rowMeanValue, positiveDetections, positiveDetectionsPerc)
  }
}

#2 nc2
negTotal2 <- data.frame(matrix(nrow=length(negZOTUs2),ncol=5))
rownames(negTotal2) <- negZOTUs2
colnames(negTotal2) <- c("readsInNeg","%ofTotalRead", "meanReads", "PresentInPositive", "%ofPositivepresent")

for (x in negZOTUs2) {
  if (x %in% negZOTUs2) {
    negValue <- freqTable[x, 'nc2']
    negValuePerc <- negValue / rowSums(freqTable[x,nc2Rows]) * 100
    rowMeanValue <- rowMeans(freqTable[x,nc2Rows])
    positiveDetections <- sum(freqTable[x, nc2Rows] > 0)-1
    positiveDetectionsPerc <- (sum(freqTable[x,nc2Rows] > 0)-1) / (length(nc2Rows)-1) * 100
    negTotal2[x,] <- c(negValue, negValuePerc, rowMeanValue, positiveDetections, positiveDetectionsPerc)
  }
}

#3 nc3
negTotal3 <- data.frame(matrix(nrow=length(negZOTUs3),ncol=5))
rownames(negTotal3) <- negZOTUs3
colnames(negTotal3) <- c("readsInNeg","%ofTotalRead", "meanReads", "PresentInPositive", "%ofPositivepresent")

for (x in negZOTUs3) {
  if (x %in% negZOTUs3) {
    negValue <- freqTable[x, 'nc3']
    negValuePerc <- negValue / rowSums(freqTable[x,nc3Rows]) * 100
    rowMeanValue <- rowMeans(freqTable[x,nc3Rows])
    positiveDetections <- sum(freqTable[x, nc3Rows] > 0)-1
    positiveDetectionsPerc <- (sum(freqTable[x,nc3Rows] > 0)-1) / (length(nc3Rows)-1) * 100
    negTotal3[x,] <- c(negValue, negValuePerc, rowMeanValue, positiveDetections, positiveDetectionsPerc)
  }
}

#4 nc4
negTotal4 <- data.frame(matrix(nrow=length(negZOTUs4),ncol=5))
rownames(negTotal4) <- negZOTUs4
colnames(negTotal4) <- c("readsInNeg","%ofTotalRead", "meanReads", "PresentInPositive", "%ofPositivepresent")

for (x in negZOTUs4) {
  if (x %in% negZOTUs4) {
    negValue <- freqTable[x, 'nc4']
    negValuePerc <- negValue / rowSums(freqTable[x,nc4Rows]) * 100
    rowMeanValue <- rowMeans(freqTable[x,nc4Rows])
    positiveDetections <- sum(freqTable[x, nc4Rows] > 0)-1
    positiveDetectionsPerc <- (sum(freqTable[x,nc4Rows] > 0)-1) / (length(nc4Rows)-1) * 100
    negTotal4[x,] <- c(negValue, negValuePerc, rowMeanValue, positiveDetections, positiveDetectionsPerc)
  }
}


# combine all 4 together
negTotal <- rbind(negTotal1,negTotal2,negTotal3,negTotal4)




# relax filtering -> for any nc sequence that is found in the sample, if the abundance is lower than x10 of the nc sequence in the nc sample -> remove. 

# 1. create a new df 
# 2. forloop to set all values less than 10 x the negative control value to 0

#nc1
freqTable.relaxedFilter <- freqTable
for (negZOTU in negZOTUs1) {
  condition <- freqTable[negZOTU, 'nc1'] * 10
  freqTable.relaxedFilter[negZOTU, nc1Rows] <- ifelse(freqTable[negZOTU, nc1Rows] > condition, 
                                                 freqTable[negZOTU, nc1Rows], 0)
}

#nc2
for (negZOTU in negZOTUs2) {
  condition <- freqTable[negZOTU, 'nc2'] * 10
  freqTable.relaxedFilter[negZOTU, nc2Rows] <- ifelse(freqTable[negZOTU, nc2Rows] > condition, 
                                                 freqTable[negZOTU, nc2Rows], 0)
}

#nc3
for (negZOTU in negZOTUs3) {
  condition <- freqTable[negZOTU, 'nc3'] * 10
  freqTable.relaxedFilter[negZOTU, nc3Rows] <- ifelse(freqTable[negZOTU, nc3Rows] > condition, 
                                                      freqTable[negZOTU, nc3Rows], 0)
}

#nc4
for (negZOTU in negZOTUs4) {
  condition <- freqTable[negZOTU, 'nc4'] * 10
  freqTable.relaxedFilter[negZOTU, nc4Rows] <- ifelse(freqTable[negZOTU, nc4Rows] > condition, 
                                                      freqTable[negZOTU, nc4Rows], 0)
}


sum(freqTable[,])
sum(freqTable.relaxedFilter[,])
sum(freqTable.relaxedFilter[,]) / sum(freqTable[,]) * 100
rownames(freqTable.relaxedFilter)[which(rowSums(freqTable.relaxedFilter[,]) == 0)]

# clean up the freq.relaxfilted table and export 
# remove negative control 
freqTable.relaxedFilter.nNeg <- freqTable.relaxedFilter[,1:108]

# remove all contaminant 
freqTable.relaxedFilter.contemRM <- 
  freqTable.relaxedFilter.nNeg[rowSums(freqTable.relaxedFilter.nNeg[,])>0,]

# write it out & convert it back to biom
write.csv(freqTable.relaxedFilter.contemRM, file = "decontaminate/freqTable.relax.contemRM.csv")
biom_obj <- make_biom(freqTable.relaxedFilter.contemRM)
write_biom(biom_obj, "decontaminate/table.relaxfilter.biom")



# clean up the rep-seq table and export 
# remove empty reps 
sequenceTableConRmv <- sequenceTable[rownames(freqTable.relaxedFilter.contemRM),]

# write it out as a new fasta file
filepath <- "E:/OneDrive - connect.hku.hk/#2 PhD/#2 research project/#1 biodiveristy/#2 CRF EXP1/#2 dataAnalysis_all/1_data/decontaminate/rep-seqs.RelaxedFilter.fasta"
writeXStringSet(sequenceTableConRmv, file = filepath, format = "fasta")
