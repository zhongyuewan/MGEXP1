# Key challenge 
# the way MACSE operates is really inpredictable, if the AA sequence alighed to is too big it wont work; if the alighnment sequence is too big it wont work.
# so in order for MACSE to work, we need to break down our sequnece into smaller files that only has 10 reads/file and then break down the alighed sequence into file <=2,500 reads
# aliged each of the 10reads seqeunce file to all of the 2,500 read alighed file and then merge all the aligned reads. 

------------------------------------------------------------------------------------------------------------------------------------------------------------
# 1/4 data preperation break down our sequence table into 10 reads/file 
## load lib
library(Biostrings)
library(dplyr)
library(biomformat)
## set WD 
setwd("filepath")
rm(list=ls())

## get data
# convert biom to csv 
sequenceTable <- readDNAStringSet('rep-seqs.RelaxedFilter.fasta')
# how long is it? 
length(sequenceTable) #36259 total reads 

# make a function to break it down and write it out 
# n sequence/file 


breakDown <- function (seqFile,n){
  cleanRowN <- floor(length(seqFile)/n)
  lastRowN <- cleanRowN +1
  for (i in 1:cleanRowN){
    
    seqTemp <- seqFile[((i-1)*n+1):(n*i)]
    
    out_file <- paste0("repSeqBy10", "_", i, ".fasta")
    writeXStringSet(seqTemp, file = out_file, format = "fasta", width = 400)
    
  }
  
  seqTemp <- seqFile[(cleanRowN*n+1):length(seqFile)]
  out_file <-  paste0("repSeqBy10", "_", lastRowN, ".fasta")
  writeXStringSet(seqTemp, file = out_file, format = "fasta", width = 400)
  
} 

breakDown(sequenceTable, 10) # after breakdown, there are 3626 .fasta files with 10 reads/file. 


------------------------------------------------------------------------------------------------------------------------------------------------------------
# 2/4 align 
# before running this, read the following so you know what's going on 
# This code is NOT THE BEST way to run one command in 6 files
# so to align all the sequence files (in the following case 3626) with a section of the aligned AA db 
# because if it's the full DB then it does not work out for some reason 
# so we cut this DB into 6 parts and put 1/6 of them in each file )
# with the software macse.jar
# so all 1-6 files will have the following:
# 1) all the rep-seq that are broken into 10 reads/file
# 2) macse.jar
# 3) mergedAll1.fasta (1-6)

## cleanup 
rm(list=ls())

# Loop through the file numbers
for (j in 1:6){# for 1~6 aligned sequence database 
  wd <- paste0("workingpath",j,"/")
 
   setwd(wd)
  
for (i in 1:3626) {
  # Construct the file names
  align_file <- paste0("mergedAll", j ,".fasta")
  seq_file <- paste0("repSeqBy10_", i, ".fasta")
  
  # Construct the command to execute
  command <- paste0(
    "java -jar macse.jar -prog enrichAlignment -align ", align_file,
    " -seq ", seq_file, " -gc_def 5 -maxSTOP_inSeq 0 -maxDEL_inSeq 3 -maxFS_inSeq 0 -maxINS_inSeq 0 -fixed_alignment_ON -max_NT_trimmed 15"
  ) # 0 frame shift, max 3 deletion, 0 insersion
  
  # Execute the command
  system(command, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  NTname <- paste0("mergedAll", j ,"_NT.fasta")
  AAname <- paste0("mergedAll", j ,"_AA.fasta")  
  Staname <- paste0("mergedAll", j ,"_stats.csv")
  
  file.rename(NTname, paste0("alignNT",i,".fasta"))
  file.rename(AAname, paste0("alignAA",i,".fasta"))
  file.rename(Staname, paste0("alignStats",i,".csv"))
}
  
}

------------------------------------------------------------------------------------------------------------------------------------------------------------
# 3/4 Retrieve aligned sequence name and remove n-aligned reads 
# after the above step 2, there will be 3626 .csv file (1 from each mini alignment) in each of the 6 AA alignment file. 
# now i just need to read all the 3626 x 6 .csv file and retrieve all the successful alignment id for downstream analysis 

## load lib
rm(list=ls())

## set WD 
for (j in 1:6){# for 1~6 aligned sequence database 
  wd <- paste0("/filepath",j,"/")
  setwd(wd)
  
# read file 
csv_files <- list.files(pattern = "*.csv", full.names = TRUE)

merged_data <- data.frame()

# Loop over each CSV file
for (file in csv_files) {
  # Read the CSV file
  data <- read.csv(file, header = TRUE)
  # Append the data to the merged_data data frame
  merged_data <- rbind(data, merged_data)  
}

#break down the big file and write into ndf
ndf <- data.frame(matrix(0, nrow=nrow(merged_data), ncol=8))
cloumn <- unlist(strsplit(names(merged_data), "\\."))
names(ndf) <- cloumn
for (i in 1:nrow(merged_data)){
  ndf[i,] <- unlist(strsplit(merged_data[i,],";"))
 
}
alignedSeq <- ndf[(ndf[,2]=="yes"),1]
length(alignedSeq)
writeLines(alignedSeq, file.path(dirname(getwd()),paste0("alignedSeq",j,".txt")))
}
------------------------------------------------------------------------------------------------------------------------------------------------------------
# 4/4 lastly move all the alignedSeqx.txt (1~6) file to the same folder and runs some quick codes to mergy all of them 


## load lib
rm(list=ls())

## set WD 
setwd("filepath")

# read file 
seq1 <- readLines("alignedSeq1.txt")
seq2 <- readLines("alignedSeq2.txt")
seq3 <- readLines("alignedSeq3.txt")
seq4 <- readLines("alignedSeq4.txt")
seq5 <- readLines("alignedSeq5.txt")
seq6 <- readLines("alignedSeq6.txt")


# merge the list 
shared <- intersect(seq1, seq2)
for(i in 3:6){
  fileName <- paste0("seq",i)
  shared <- intersect(shared, get(fileName))
}

uniSeq <- union(seq1, seq2)
for(i in 3:6){
  fileName <- paste0("seq",i)
  uniSeq <- union(uniSeq, get(fileName))
}


# write it out 
writeLines(uniSeq, "cleanSeq.txt") 
# Collect the sequence that got picked up by any of the 1~6 aligned files. 21,757 reads aligned to the database. 



