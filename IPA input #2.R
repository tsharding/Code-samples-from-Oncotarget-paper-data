### for IPA input 10_16_17

#input all diff exp (exclude temoral comparisons)
#add input for unique to combo
#Filter at |FC|> 1.5 and two -> export to separate folders for input to IPA

IPAinput2.names.MMM1 <- c("MA4vME4","MA5vME5","MA5vMB5","MA5vMG5","MMM1_UTC")
IPAinput2.names.FLAM76 <- gsub("M","F",IPAinput2.names.MMM1)
IPAinput2.names.FLAM76 <- gsub("FFF1","FLAM76",IPAinput2.names.FLAM76)

UTC.func <-function(line){ #accepts "M" or "F" - returns xA5vxG5 of filtered for UTC at |fc|>1.5 for the respective line
  #get
  E <- get(gsub("x",line,"xA5vxE5"))
  B <- get(gsub("x",line,"xA5vxB5"))
  G <- get(gsub("x",line,"xA5vxG5"))
  #filter
  EE <- filter.df.data(E,1.5)[,"gene"]
  BB <- filter.df.data(B,1.5)[,"gene"]
  GG <- filter.df.data(G,1.5)[,"gene"]
  #get UTC genes
  UTC.genes <- setdiff(GG,union(EE,BB))
  print(length(UTC.genes))
  return(G[which(G$gene %in% UTC.genes),])
  }


MMM1_UTC <- UTC.func("M")
FLAM76_UTC <- UTC.func("F")

#loop to create files for import #2
setwd("~/Data after 11_11_15/big exp/RNA-seq/IPA Import 2")
for(x in c(IPAinput2.names.MMM1,IPAinput2.names.FLAM76)){
  print(x)
  DF <- get(x)
  write.table(filter.df.data(DF,1.5), file = paste(x,"1.5.txt",sep = "_"), sep = '\t',quote = FALSE,eol = "\r", row.names = FALSE)
  write.table(filter.df.data(DF,2), file = paste(x,"2.txt",sep = "_"), sep = '\t',quote = FALSE,eol = "\r", row.names = FALSE)
}
rm(x,DF)