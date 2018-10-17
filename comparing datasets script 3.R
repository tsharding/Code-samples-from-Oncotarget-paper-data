#generate table of genes most enriched in all data sets
#exclude genes not identified in my data

#get gene names
source("http://bioconductor.org/biocLite.R")
biocLite("AnnotationDbi")
biocLite("org.Hs.eg.db")
library("AnnotationDbi")#use this package to convert entrez IDs to gen symbols
library("org.Hs.eg.db")
IDs.mydata.GENENAMES <- mapIds(org.Hs.eg.db, keys=IDs.mydata, column="GENENAME", keytype="SYMBOL", multiVals="first")
#remove commas from gene names so data can be exported to .csv
IDs.mydata.GENENAMES <- unlist(lapply(IDs.mydata.GENENAMES,gsub,pattern = ",",replacement = ";"))

#initiate data.frame with rownames as all RNAseq IDs
My.Data.DF <- data.frame("GENE NAME" = IDs.mydata.GENENAMES,row.names = IDs.mydata)
My.Data.DF$total.freq <- 0
My.Data.DF$total.freq.UP <- 0
My.Data.DF$total.freq.DOWN <- 0



#lists of dataset names for certain colums
FLAM_EPZ <- c("FA1vFE1","FA4vFE4","FA5vFE5")
My.Data.DF$FLAM_EPZ <- 0
MMM1_EPZ <- c("MA1vME1","MA4vME4","MA5vME5")
My.Data.DF$MMM1_EPZ <- 0
BOTH_EPZ <- c(FLAM_EPZ,MMM1_EPZ)
My.Data.DF$BOTH_EPZ <- 0

for(x in MyData_names){
  print(x,quote =  FALSE)
  UP <- get(paste(x,"UP",sep = "_"))
  DOWN <- get(paste(x,"DOWN",sep = "_"))
  BOTH <- c(UP,DOWN)
  #all of my data
  My.Data.DF$total.freq <- My.Data.DF$total.freq + (row.names(My.Data.DF) %in% BOTH)
  My.Data.DF$total.freq.UP <- My.Data.DF$total.freq.UP + (row.names(My.Data.DF) %in% UP)
  My.Data.DF$total.freq.DOWN <- My.Data.DF$total.freq.DOWN + (row.names(My.Data.DF) %in% DOWN)
  #Both epz
  if(x %in% BOTH_EPZ){
    My.Data.DF$BOTH_EPZ <- My.Data.DF$BOTH_EPZ + (row.names(My.Data.DF) %in% BOTH)
  }
  #FLAM EPZ
  if(x %in% FLAM_EPZ){
    My.Data.DF$FLAM_EPZ <- My.Data.DF$FLAM_EPZ + (row.names(My.Data.DF) %in% BOTH)
  }
  #MMM1 EPZ
  if(x %in% MMM1_EPZ){
    My.Data.DF$MMM1_EPZ <- My.Data.DF$MMM1_EPZ + (row.names(My.Data.DF) %in% BOTH)
  }
}
My.Data.DF.sorted <- My.Data.DF[order(-My.Data.DF$total.freq),]
My.Data.DF.sorted <- My.Data.DF.sorted[My.Data.DF.sorted$total.freq != 0,]


write.csv(My.Data.DF.sorted,"My data reccouring genes sorted.csv",quote = FALSE)


##################################################################################################################
#enrichment of other datasets in my data
Enrichment.negative.DF  <- My.Data.DF[My.Data.DF$total.freq == 0,1:2] #not sigdiff in my data ever
Enrichment.negative.DF$total.freq.UP <- 0
Enrichment.negative.DF$total.freq.DOWN <- 0

Enrichment.DF <- My.Data.DF[My.Data.DF$total.freq != 0,1:2] #sigdiff in my data at least once
Enrichment.DF$total.freq.UP <- 0
Enrichment.DF$total.freq.DOWN <- 0

#includes genes shared between sisters for respective conditions
Enrichment.E.common.DF <- My.Data.DF[FA5vFE5_common_with_sister_lowFC,1:4]
Enrichment.G.common.DF <- My.Data.DF[FA5vFG5_common_with_sister_lowFC,1:4]
Enrichment.BtoG.common.DF <- My.Data.DF[FB5vFG5_common_with_sister_lowFC,1:4]
Enrichment.EtoG.common.DF <- My.Data.DF[FE5vFG5_common_with_sister_lowFC,1:4]

#clean/initiate above df's
Enrichment.names <- c("Enrichment.negative.DF","Enrichment.DF","Enrichment.E.common.DF","Enrichment.G.common.DF","Enrichment.BtoG.common.DF","Enrichment.EtoG.common.DF")
for(x in Enrichment.names){
  print(x,quote=FALSE)
  y<-get(x)
  y$total.freq.NEITHER <- 0
  y$MyDataRecourance <-y$total.freq
  y$total.freq.Datasets <- ""
  y$total.freq <- 0
  y$total.freq.UP <- 0
  y$total.freq.DOWN <- 0
  assign(x,y)
}
rm(x,y)


#set lists of datasets to be compared to; checked to be in ls(), example: >sum(sapply(Comparison.names.DOWN,function(x) length(get(x))) == 0)
Bayer.list.up <- unlist(lapply(Bayer.list,function(x) paste("Bayer",x,"up",sep = ".")))
Bayer.list.down <- unlist(lapply(Bayer.list,function(x) paste("Bayer",x,"down",sep = ".")))

Comparison.names.UP <- c(Bayer.list.up,"MONOMETH.up",oncomine.names[1],unlist(lapply(MCCABE.filenames,paste,"up",sep = ".")),"INA6.up","IDs.ma2.INA6up")
#remove empty lists
Comparison.names.UP <- Comparison.names.UP[sapply(Comparison.names.UP,function(x) length(get(x))) != 0]

Comparison.names.DOWN <- c(Bayer.list.down,"MONOMETH.down",unlist(lapply(MCCABE.filenames,paste,"down",sep = ".")),"INA6.down","IDs.ma2.INA6down")
#remove empty lists
Comparison.names.DOWN <- Comparison.names.DOWN[sapply(Comparison.names.DOWN,function(x) length(get(x))) != 0]

Comparison.names.NEITHER <-c(oncomine.names[c(2,6:10)],"MM_unique_H3K27me3_targets","MM_unique_bivalent_targets")

#master loop for running comparisons #######################################################
for(x in Enrichment.names){
  print(x, quote = FALSE)
  print("",quote = FALSE)
  x.names <- row.names(get(x))
  x.data <- get(x)
  for(y in Comparison.names.UP){
    print(y, quote = FALSE)
    z <- (x.names %in% get(y))
    x.data$total.freq <- x.data$total.freq + z
    x.data$total.freq.UP <- x.data$total.freq.UP + z
    x.data$total.freq.Datasets[z] <- paste(x.data$total.freq.Datasets[z],y,sep = " ; ")
  }
  print("",quote = FALSE)
  for(y in Comparison.names.DOWN){
    print(y, quote = FALSE)
    z <- (x.names %in% get(y))
    x.data$total.freq <- x.data$total.freq + z
    x.data$total.freq.DOWN <- x.data$total.freq.DOWN + z
    x.data$total.freq.Datasets[z] <- paste(x.data$total.freq.Datasets[z],y,sep = " ; ")
  }
  print("",quote = FALSE)
  for(y in Comparison.names.NEITHER){
    print(y, quote = FALSE)
    z <- (x.names %in% get(y))
    x.data$total.freq <- x.data$total.freq + z
    x.data$total.freq.NEITHER <- x.data$total.freq.NEITHER + z
    x.data$total.freq.Datasets[z] <- paste(x.data$total.freq.Datasets[z],y,sep = " ; ")
  }
  assign(x,x.data)
  print("",quote = FALSE)
}
rm(x,y,z,x.names,x.data)

#sort each df and save to csv
setwd("~/Data after 11_11_15/EZH2 informatics from literature/enrichment csv files - comparing to my data")
Enrichment.file.name.keys <- c("negative.csv","all data.csv","A5vE5.csv","A5vG5.csv","B5vG5.csv","E5vG5.csv")
names(Enrichment.file.name.keys)<-Enrichment.names
for(x in Enrichment.names){
  x.data <- get(x)
  #sort by total.freq
  x.data <- x.data[order(-x.data$total.freq),]
  assign(paste(x,"sorted", sep = "."),x.data)
  #write as .csv
  write.csv(x.data,Enrichment.file.name.keys[x],quote = FALSE)
}
rm(x,x.data)



##################################################
#for comparing datasets to Morgan 2017 paper
Morgan.2017.up.KMS11 <- read.table(file = "clipboard",sep = "\t",stringsAsFactors = FALSE, quote = "")[,1]
Morgan.2017.up.KMM1 <- read.table(file = "clipboard",sep = "\t",stringsAsFactors = FALSE, quote = "")[,1]
Morgan.2017.down.KMS11 <- read.table(file = "clipboard",sep = "\t",stringsAsFactors = FALSE, quote = "")[,1]
Morgan.2017.down.KMM1 <- read.table(file = "clipboard",sep = "\t",stringsAsFactors = FALSE, quote = "")[,1]
Morgan.2017.list <- c("Morgan.2017.up.KMM1","Morgan.2017.up.KMS11","Morgan.2017.down.KMS11","Morgan.2017.down.KMM1")

Enrichment.sorted.names.list <- sapply(Enrichment.names,paste,"sorted",sep = ".")
  
setwd("~/Data after 11_11_15/EZH2 informatics from literature/enrichment csv files - comparing to my data")  
for(x in Morgan.2017.list){
  print(x,quote = FALSE)
  y <- get(x)
  x.data <- get(Enrichment.sorted.names.list[2]) #gets sorted df for all data
  x.data <- x.data[row.names(x.data) %in% y,]
  print(length(y),quote = FALSE)
  print(sum(row.names(x.data) %in% y),quote = FALSE)
  #write as .csv
  write.csv(x.data,paste(x,"csv",sep = "."),quote = FALSE)
}
rm(y,x,x.data)

##########################
#subset enriched data for IPA import


