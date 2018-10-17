####finding good housekeeping gene

MyData_names

HKG <- data.frame(GENENAME = My.Data.DF[,1],KEEP = 0,MEAN.FPKM.MMM1 = 0,MEAN.FPKM.FLAM = 0,MAX.FC=0,row.names = row.names(My.Data.DF))
# 0 = keep 

for(x in MyData_names){
  y <- get(x)
  HKG$KEEP <- HKG$KEEP + (y$q_value < 0.05)
  HKG$KEEP <- HKG$KEEP + (y$value_1 < 1)
  HKG$KEEP <- HKG$KEEP + (y$value_2 < 1)
  #distinguish MMM1 and FlAM conditions within loop
  if (grepl("M",x) == TRUE) {HKG$MEAN.FPKM.MMM1 <- HKG$MEAN.FPKM.MMM1 + y$value_1} #use value_1 value only
  if (grepl("F",x) == TRUE) {HKG$MEAN.FPKM.FLAM <- HKG$MEAN.FPKM.FLAM + y$value_1} #use value_1 value only
  MAX <- (abs(y$log2.fold_change.)>abs(HKG$MAX.FC))
  HKG$MAX.FC[MAX] <- y$log2.fold_change.[MAX]
}
rm(x,y)
HKG$MEAN.FPKM.MMM1 <- HKG$MEAN.FPKM.MMM1/11
HKG$MEAN.FPKM.FLAM <- HKG$MEAN.FPKM.FLAM/11
HKG <- HKG[HKG$KEEP == 0,]
HKG <- HKG[order(abs(HKG$MAX.FC)),]

#write to .csv
setwd("~/Data after 11_11_15/EZH2 informatics from literature")
write.csv(HKG,"Housekeeping genes based on my data.csv",quote = FALSE)
