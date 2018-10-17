# Analysis of microarray data from paper:
# Run with R v3.2.4
# "EZH2 Inhibition Blocks Multiple Myeloma Cell Growth through Upregulation of Epithelial Tumor Suppressor Genes"
#   (2015)Hernando et al - Bayer Pharm
#   DOI: 10.1158/1535-7163.MCT-15-0486 
#*****authors used 'ArrayExpress' package 

# data downloaded from http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3540/
#     3 zip files: 	E-MTAB-3540.raw.1.zip, E-MTAB-3540.raw.2.zip, E-MTAB-3540.raw.3.zip
#           - 95 '.CEL' files
#           - array: pd.hugene.2.1.st (plate info: http://www.affymetrix.com/catalog/131554/AFFY/Human+Gene+ST+Array+Plates#1_1)
#     combined in folder wd:"~/Data after 11_11_15/EZH2 informatics from literature/Data from Bayer Pham paper _ Microarray EZH2i kinetic GEPs/Combined raw data"
#     some file numbers between 1 and 100 missing: 22,35,54,62,93 - presumed as missing replicates (5 per condition)

# set WD for raw data 
  setwd("~/Data after 11_11_15/EZH2 informatics from literature/Data from Bayer Pham paper _ Microarray EZH2i kinetic GEPs/Combined raw data")

# make list of raw data file names
  rawdatafilenames <- list.celfiles(full.names=FALSE)

# make phenodata data frame
  # row names are file name (Fromat: "###_#_celline_treatment.CEL")
  # sample numbers 1-100 (missing numbers: 22,35,54,62,93)
  # cell_lines: "KMS11","KMS12PE","KMS28BM","KMS34","LP1","MOLP8","NCIH929","RPMI8226","KMS20","U266B1"
  # treatments: "DMSO" or "E7438" (4 or 5 samples each)
  file_elements <- do.call(rbind,strsplit(rawdatafilenames,"_")) #splits file names by "_" and converts to matrix
  pd <- data.frame(row.names = factor(rawdatafilenames), sample = file_elements[,1], cell_line = file_elements[,3], treatment = sub(".CEL","",file_elements[,4]))

# data pre-processing via oligo package [should also be compatible with the 'xps' package]
  # used Carvalho protocol found online:
  # http://www.bioconductor.org/packages//2.7/bioc/vignettes/oligo/inst/doc/V5ExonGene.pdf
  # 'oligo' package user guide: https://www.bioconductor.org/packages/3.3/bioc/vignettes/oligo/inst/doc/oug.pdf
  source("http://bioconductor.org/biocLite.R")
  biocLite("oligo")
  library("oligo")
      # biocLite("pd.hugene.2.0.st") #'oligo' will automatically import this affymatrix annotation package
      # library("pd.hugene.2.0.st")

  # read .CEL files
  dat <- read.celfiles(list.celfiles())
  
  # pre-process data via rma() function - returns ExressionSet class data
  #       - Background correction, Normalization (meadian polish) & Expression calculation
  ppData <- rma(dat) #default is transcript level data: rma(,target = "core")
  
  # assign pheno data columns to ExpressionSet object
  ppData[['cell_line']] <- pd$cell_line
  ppData[['sample']] <- pd$sample
  ppData[['treatment']] <- pd$treatment
  
  #annotate ExpressionSet object from rma() output with gene-symbol info as featurdata(eset)
  #note: getNetAffx() from oligo package does not work here
  source("http://bioconductor.org/biocLite.R")
  biocLite("annotate")
  library("annotate")
  biocLite("hugene21sttranscriptcluster.db") #don't use "pd.hugene.2.0.st" file or "HuGene-2_1-st-v1.na36.hg19.transcript" file for annotation
  library("hugene21sttranscriptcluster.db")
  #add gene-symbol info to featuredata(eset)
  #examples: <https://support.bioconductor.org/p/63834/>
  #         <http://www.gettinggeneticsdone.com/2012/01/annotating-limma-results-with-gene.html>
  #note: IPA does not need annotation to mach transcript IDs to symbol
    ID <- featureNames(ppData)
    Symbol <- getSYMBOL(ID,"hugene21sttranscriptcluster.db")
    fData(ppData) <- data.frame(Symbol=Symbol)
    #24082 transcripts are NA
    #29535 transcripts assigned symbol
  
  
   
# Create expression matrix to check expression distribution
  expression_matrix <- exprs(ppData)
  # plot histogram of mean expression for each 
  hist(rowMeans(expression_matrix),breaks=1000, main="Figure 1: Histogram of expression means accross all samples",xlab="Normalized expression", ylab="Number of probes")

# Differential expression analysis
  # Bayer paper used paired t-test with FDR (Benjamini-Hochberg)<0.1 and filter by fold change >1.5 using Expressionist-GeneData software
  
  #loads genefilter library which contains ttest function
  #NOT USED FOR ANALYSIS - CONSIDER FOR FILTERING FUNCTIONS
#      biocLite("genefilter")
#      library("genefilter")
#      #ttest based on 'treatment' phenodata in KMS11
#      tt <- rowttests(ppData[,ppData$cell_line == 'MOLP8'], "treatment")
#      #p-value ajustment
#      pBH <- p.adjust(tt$p.value, method="BH")
#      tt_BHfilter <- (tt[which(pBH<0.01),])
  

#subset ppData for testing in one line (KMS11)  
  cell.line <- ppData[,ppData$cell_line=="KMS11"]
  
#limma differential expression  
  #LIMMA USER GUIDE: <http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf>
  source("http://bioconductor.org/biocLite.R")
  biocLite("limma")
  library("limma")
  
  #create function that takes cell line (one from pd$cell_line) and returns results data.frame
  expression.analysis <- function(x){
    # x = name of cell line from pd$cell_line for subsetting ppData 
      require(limma)
      cell.line <- ppData[,ppData$cell_line == x]
      print(paste("cell line:",x))
    #differential expression via limma; modern t-test
      design <- model.matrix(~factor(cell.line$treatment)) #key for differential expression: 0=DMSO 1=E6438
      fit <- lmFit(cell.line, design)
      ebayes <- eBayes(fit) #modern t-test
      assign(paste(x,"_diffexprs_ebayes", sep = ""),ebayes, envir=.GlobalEnv)
    #calculate mean expression for each condition per transcript
      CL.exprs <- exprs(cell.line)
      DMSO.mean <- rowMeans(CL.exprs[,cell.line$treatment == "DMSO"])
      E6438.mean <- rowMeans(CL.exprs[,cell.line$treatment == "E7438"])
    #create new data frame with all results and info
      #make ebayes data.frame for easy subsetting, also perform "BH" p.value adjustment
      df.ebayes <- as.data.frame(ebayes)
      
      #fix colnames
      cndf <- colnames(df.ebayes)
      cndf[2] <- "logFC"
      #note logFC is calculated here based on geometric mean not arithmetic mean: <https://support.bioconductor.org/p/42922/>
      cndf[8] <- "Mean expression"
      cndf[14] <- "p.value"
      colnames(df.ebayes) <- cndf
      
      #multiple testing correction : "BH" 
      df.ebayes$q.value <- p.adjust(df.ebayes$p.value, method = "BH")
      
      #combine relavent data into data.frame
      df <- data.frame(cbind(fData(cell.line),DMSO.mean,E6438.mean,df.ebayes$`Mean expression`,df.ebayes$logFC,df.ebayes$p.value,df.ebayes$q.value))
      
      #fix colnames again
      cndf <- colnames(df)
      cndf[4]<-"total.mean"
      cndf[5]<-"logFC"
      cndf[6]<-"p.value"
      cndf[7]<-"q.value"
      colnames(df) <- cndf
      
      assign(paste(x,"_diffexprs_results", sep = ""),df, envir=.GlobalEnv)
      
      return(df)
      
    rm(x,cell.line,design,fit,ebayes,CL.exprs,DMSO.mean,E6438.mean,df,df.ebayes,cndf)
    }
  
  #write diffxprs data to tab delimited text file for import into IPA, accepts data.frame - writes to .txt file
  write.results <- function(y,CL.name){
      write.table(y, file = paste(CL.name,"_results.txt",sep=""), sep = '\t',quote = FALSE,eol = "\r", col.names = NA,row.names = TRUE)
      rm(y)
  }
  
  #loops through all cell lines. For each: diff expression analysis w/ "BH" p.value adjustment, assigns 'xxx_diffexprs_ebayes' and 'xxx_diffexprs_results', returns data.frame of results, writes data.frame to .txt file
  for (z in unique(pd$cell_line)) {
    results <- expression.analysis(z)
    write.results(results,z)
    rm(results)
  }
  

  