# every ARMS contains organic materials of 3 differet fractions -> 500um/106um/sessile
# the follow code combainds all of the three fractions by ARMS, producing the following new files: 
# metadata by ARMS & frequency table by ARMS

## set WD 
# mac
setwd("filepath")
rm(list=ls())


## load lib
library(dplyr)
library(ggplot2)
library(Biostrings)
library(biomformat)


## get data
metaData <- read.table('sample-metadata.tsv', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, comment.char = '')
freqTable <- read.csv('decontaminate/clean/cluster/feaTable97.csv', header = TRUE, row.names = 1)
sequenceTable <- readDNAStringSet('decontaminate/clean/cluster/dna-sequences97.fasta')


# gsub X
names(freqTable) <- gsub("X","",names(freqTable))

## remove the 3 duplicate data point and the negative condtol
sampleRow <- row.names(metaData[metaData[,2]=='f' & metaData[,1] != 'na',])
#check how long is the list, should be 105 because of 35 arms with 3 fractions each
length(sampleRow) == 105 

## make a new metaData & new frequency table, no need to make new sequence table
# new metaData d1 <- trim down wiht the sampleRow
metaDataByArms.d1 <- metaData[sampleRow,]

# new metaData d2 <- combine the read counts of all three fractions 
add_summed_columns <- function(df) {
  # Calculate the sums for columns 26, 27, 29, 30, and 32  grouped by column 1
  sums <- aggregate(df[, c(26, 27, 29, 30, 32)], by = list(df[, 1]), FUN = sum, na.rm = TRUE)
  colnames(sums) <- c("ARMSno", "sum_input", "sum_filtered", "sum_denoised", "sum_merged", "sum_nChimeric")
  
  # Merge the sums back into the original data frame 
  # merge by by.x and by.y, in the following example, both of them are armns numbers. 
  # all.x = TRUE makes it write on all the original data frame
  df <- merge(df, sums, by.x = colnames(df)[1], by.y = "ARMSno", all.x = T)
  
  return(df)
}

metaDataByArms.d2 <- add_summed_columns(metaDataByArms.d1)
metaDataByArms.d2$nChimericByInput <-metaDataByArms.d2[,38]/metaDataByArms.d2[,34] 


# now we can just take the sessile (or 106/500) out as the final dataSet
metaDataByArms.d3 <- metaDataByArms.d2[metaDataByArms.d2[,5]=="sessile",]
# clean it up
metaDataByArms <- metaDataByArms.d3[,c(1,6:22,34:39)] # because "sessile" reads are the same with "ARMS" reads

# combining a new frequency table
freqTable.d1 <- freqTable

# remove all the _2_106 data (total 3 of them)
removeColumn <- function(dataframe) {
  columns_to_remove <- grep("_2_", colnames(dataframe))
  dataframe <- dataframe[, -columns_to_remove]
  
  return(dataframe)
}

freqTable.d2 <- removeColumn(freqTable.d1)
ncol(freqTable.d2)==105 ## after clean up it should have 105 column 

# add all the rows together by arms 

groupColumns <- function(dataframe) {
  groups <- unique(substr(names(dataframe), 1, 2))
  grouped_data <- data.frame(matrix(nrow =nrow(dataframe)))
  row.names(grouped_data) <- row.names(dataframe)
  
  for (group in groups) {
    group_columns <- grep(paste0("^", group), names(dataframe), value = TRUE)
    grouped_data[group] <- rowSums(dataframe[, group_columns])
  }
  
  return(grouped_data[-1])
}

freqTableByArms <- groupColumns(freqTable.d2)



## write them out
# metadata
names(metaDataByArms)[1] <- 'featureid'
write.table(metaDataByArms, file = "decontamByArms/sample-metadata.byArms.tsv", row.names=F)


# frequence table 
# write it out & convert it back to biom
write.csv(freqTableByArms, file = "decontamByArms/freqTableCleanByArms.csv")
biom_obj <- make_biom(freqTableByArms)
write_biom(biom_obj, "decontamByArms/table.clean.byArms.biom")


