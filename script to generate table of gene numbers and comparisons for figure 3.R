# generat shared # of genes in sister data set for Pie charts (figure 3)

Flam <- unlist(FLAM_files[c(2,4,5:7)])
MMM1 <- unlist(MMM1_files[c(2,4,5:7)])
shared_names <- unlist(gsub("F","S",Flam))
time_vs_names <- c("F1v4","F4v5","M1v4","M4v5","S1v4","S4v5")

DF.gene.numbers <- data.frame(row.names = c(Flam,shared_names,MMM1,time_vs_names))
DF.gene.numbers$up1.5 <- 0
DF.gene.numbers$down1.5 <- 0
DF.gene.numbers$up2 <- 0
DF.gene.numbers$down2 <- 0

filter.df.data <- function(DF,FC_threshold){ #filters data for q<0.05, |FC| cuttof >= argument and FPKM >=1 in at least 1 of two values (CTRL vs DRUG)
  return(subset(DF, q_value < 0.05 & abs(log2.fold_change.) >= log2(FC_threshold) & (value_1 >= 1 | value_2 >= 1)))
}



Gene.num.data.crunch <- function(DF){ #function to fill DF.gene.numbers data frame
  for(x in c(1:5)){
    
    F<-get(Flam[x])
    M<-get(MMM1[x])
    
    F_1.5 <- filter.df.data(F,1.5)
    F_2 <- filter.df.data(F,2)
    M_1.5 <- filter.df.data(M,1.5)
    M_2 <- filter.df.data(M,2)

    F_1.5_up <- subset(F_1.5 , log2.fold_change.>0)[,"gene"]
    F_2_up <- subset(F_2 , log2.fold_change.>0)[,"gene"]
    M_1.5_up <- subset(M_1.5 , log2.fold_change.>0)[,"gene"]
    M_2_up <- subset(M_2 , log2.fold_change.>0)[,"gene"]
    F_1.5_down <- subset(F_1.5 , log2.fold_change.<0)[,"gene"]
    F_2_down <- subset(F_2 , log2.fold_change.<0)[,"gene"]
    M_1.5_down <- subset(M_1.5 , log2.fold_change.<0)[,"gene"]
    M_2_down <- subset(M_2 , log2.fold_change.<0)[,"gene"]
    
    DF[Flam[x],1] <- length(F_1.5_up)
    DF[Flam[x],2] <- length(F_1.5_down)
    DF[Flam[x],3] <- length(F_2_up)
    DF[Flam[x],4] <- length(F_2_down)
    DF[MMM1[x],1] <- length(M_1.5_up)
    DF[MMM1[x],2] <- length(M_1.5_down)
    DF[MMM1[x],3] <- length(M_2_up)
    DF[MMM1[x],4] <- length(M_2_down)

    DF[shared_names[x],1] <- sum(F_1.5_up %in% M_1.5_up)
    DF[shared_names[x],2] <- sum(F_1.5_down %in% M_1.5_down)
    DF[shared_names[x],3] <- sum(F_2_up %in% M_2_up)
    DF[shared_names[x],4] <- sum(F_2_down %in% M_2_down)
      
    if(x == 1){
      ff_1.5_up <- F_1.5_up
      ff_2_up <- F_2_up
      mm_1.5_up <- M_1.5_up
      mm_2_up <- M_2_up
      ss_1.5_up <- F_1.5_up[F_1.5_up %in% M_1.5_up]
      ss_2_up <- F_2_up[F_1.5_down %in% M_1.5_down]
      ff_1.5_down <- F_1.5_down
      ff_2_down <- F_2_down
      mm_1.5_down <- M_1.5_down
      mm_2_down <- M_2_down
      ss_1.5_down <- F_1.5_down[F_2_up %in% M_2_up]
      ss_2_down <- F_2_down[F_2_down %in% M_2_down]
      
    }
    
    if(x == 2 | x== 4){
      if(x==2) {y<-1}
      if(x==4) {y<-2}
      #flam
      DF[time_vs_names[(y)],1] <- sum(F_1.5_up %in% ff_1.5_up)
      DF[time_vs_names[(y)],2] <- sum(F_1.5_down %in% ff_1.5_down)
      DF[time_vs_names[(y)],3] <- sum(F_2_up %in% ff_2_up)
      DF[time_vs_names[(y)],4] <- sum(F_2_down %in% ff_2_down)
      #mmm1
      DF[time_vs_names[(y+2)],1] <- sum(M_1.5_up %in% mm_1.5_up)
      DF[time_vs_names[(y+2)],2] <- sum(M_1.5_down %in% mm_1.5_down)
      DF[time_vs_names[(y+2)],3] <- sum(M_2_up %in% mm_2_up)
      DF[time_vs_names[(y+2)],4] <- sum(M_2_down %in% mm_2_down)
      #shared
      DF[time_vs_names[(y+4)],1] <- sum(ss_1.5_up %in% F_1.5_up[F_1.5_up %in% M_1.5_up])
      DF[time_vs_names[(y+4)],2] <- sum(ss_1.5_down %in% F_1.5_down[F_2_up %in% M_2_up])
      DF[time_vs_names[(y+4)],3] <- sum(ss_2_up %in% F_2_up[F_1.5_down %in% M_1.5_down])
      DF[time_vs_names[(y+4)],4] <- sum(ss_2_down %in% F_2_down[F_2_down %in% M_2_down])
    }
    
    if(x==2){
      ff_1.5_up <- F_1.5_up
      ff_2_up <- F_2_up
      mm_1.5_up <- M_1.5_up
      mm_2_up <- M_2_up
      ss_1.5_up <- F_1.5_up[F_1.5_up %in% M_1.5_up]
      ss_2_up <- F_2_up[F_1.5_down %in% M_1.5_down]
      ff_1.5_down <- F_1.5_down
      ff_2_down <- F_2_down
      mm_1.5_down <- M_1.5_down
      mm_2_down <- M_2_down
      ss_1.5_down <- F_1.5_down[F_2_up %in% M_2_up]
      ss_2_down <- F_2_down[F_2_down %in% M_2_down]
      
    }
  }
  return(DF)
}

DF.gene.numbers <- Gene.num.data.crunch(DF.gene.numbers)

setwd("~/Data after 11_11_15/EZH2 informatics from literature")
write.csv(DF.gene.numbers,"Gene numbers and comparisons for figure 3 pie charts.csv",quote = FALSE)



