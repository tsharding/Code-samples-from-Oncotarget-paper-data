#evaluating distribution of up regulated and downregulated genes. 


#do downregualtion of highly expressed genes reflect decresed representation in reads or real expression changes
#determine the sum of total FPKM change for upregulated and downregulated genes

#create function for analyze a condition; input is string for condition name
UpDown.analysis <- function(condition){
    data <- get(condition)
    print("",quote = FALSE)
    print(condition,quote = FALSE)
    print(paste("total FPKM value_1:",round(sum(data$value_1),digits = 2),"- total FPKM value_2:",round(sum(data$value_2),digits = 2),sep = " "),quote = FALSE)
    
    #filter for q<0.5 (no FC cuttoff)
    data_filtered <- subset(data, q_value <= 0.05)
    up <- subset(data_filtered,log2.fold_change.>0)
    down <- subset(data_filtered,log2.fold_change.<0)
    print(paste("# of filetered genes - total(",length(data_filtered$q),") up(",length(up$q),") down(",length(down$q),")",sep = ""),quote = FALSE)
    upDeltaFPKM <- round(sum(up$value_2 - up$value_1),2)
    downDeltaFPKM <- round(sum(down$value_1 - down$value_2),2)
    totalDeltaFPKM <- round(sum(abs(data_filtered$value_2 - data_filtered$value_1)),2)
    print(paste("delta FPKM - total(",totalDeltaFPKM,") up(",upDeltaFPKM, ") down(",downDeltaFPKM,") down%(",round(100*downDeltaFPKM/totalDeltaFPKM,2),")",sep = ""),quote = FALSE)
}

for(x in MyData_names){UpDown.analysis(x)}





#####################################################
#below is pasted from other script: "cuffdiff out put - analysis of unfiltered data"
#generating customized volcano plots to veiw disribution of expression changes
#####################################################

source("http://bioconductor.org/biocLite.R")
install.packages("ggplot2")
library(ggplot2)
#creating plots to look at distribution of expression changes based on CTRL expression



#####################################################
#altered from pasted script below
#####################################################
#create function that plots on log1p scale for delta_FPKM against actual fpkm values
#with distinction between significant an not significant
plot.actual_expression_log_sig_DF <- function(DF){
  require(ggplot2)
  DFtitle <- deparse(substitute(DF))
  DF$delta_FPKM <- log(1+abs(DF$value_2 - DF$value_1),2) * (DF$value_2 - DF$value_1)/abs(DF$value_2 - DF$value_1)
  sig_index <- which(DF$significant == "yes")
  notsig_index <- which(DF$significant == "no")
  g = ggplot() +
    geom_point(aes(DF$value_1[sig_index],DF$delta_FPKM[sig_index]),size = 1, alpha=0.2,color = "red") +
    geom_point(aes(DF$value_1[notsig_index],DF$delta_FPKM[notsig_index]),size = 1, alpha=0.2,color = "black") +
    geom_smooth(color = "red", aes(DF$value_1[sig_index],DF$delta_FPKM[sig_index])) +
    geom_smooth(color = "blue", aes(DF$value_1,DF$delta_FPKM)) +
    geom_hline(yintercept = 0,color = 'red')+
    geom_vline(xintercept = 1,color = 'green')+
    geom_vline(xintercept = 10,color = 'green')+
    geom_vline(xintercept = 100,color = 'green')+
    geom_vline(xintercept = 1000,color = 'green')+
    geom_vline(xintercept = 10000,color = 'green')+
    labs(x = "CTRL expression(FPKM)", y = "log2(1+delta_FPKM)", title = DFtitle)+
    scale_x_log10()+
    theme(legend.position = "none")+
    theme_bw()
  return(g)
}

plot.actual_expression_log_sig_DF(subset(FA5vFG5, q_value <= 0.05 & value_1 < 1000))
plot.actual_expression_log_sig_DF(subset(MA5vMG5, q_value <= 0.05))


################ change sig to include 1.5 |FC| cutoff 
plot.actual_expression_log_sig_DF_FC <- function(DF){
  require(ggplot2)
  DFtitle <- deparse(substitute(DF))
  DF$delta_FPKM <- log(1+abs(DF$value_2 - DF$value_1),2) * (DF$value_2 - DF$value_1)/abs(DF$value_2 - DF$value_1)
  sig_index <- which(DF$significant == "yes" & abs(DF$log) >= 0.584963 & DF$value_1 >= 1)
  notsig_index <- which((DF$significant == "yes" & abs(DF$log) < 0.584963)|(DF$significant == "yes" & DF$value_1 < 1))
  g = ggplot() +
    geom_point(aes(DF$value_1[sig_index],DF$delta_FPKM[sig_index]),size = 1, alpha=0.2,color = "red") +
    geom_point(aes(DF$value_1[notsig_index],DF$delta_FPKM[notsig_index]),size = 1, alpha=0.2,color = "black") +
    geom_smooth(color = "red", aes(DF$value_1[sig_index],DF$delta_FPKM[sig_index])) +
    geom_smooth(color = "blue", aes(DF$value_1[union(sig_index,notsig_index)],DF$delta_FPKM[union(sig_index,notsig_index)])) + #fixed to have q_value cutoff
    geom_hline(yintercept = 0,color = 'red')+
    geom_vline(xintercept = 1,color = 'green')+
    geom_vline(xintercept = 10,color = 'green')+
    geom_vline(xintercept = 100,color = 'green')+
    geom_vline(xintercept = 1000,color = 'green')+
    geom_vline(xintercept = 10000,color = 'green')+
    labs(x = "CTRL expression(FPKM)", y = "log2(1+delta_FPKM)", title = DFtitle)+
    scale_x_log10()+
    theme(legend.position = "none")+
    theme_bw()
  return(g)
}

#### ranked expr
plot.ranked_expression_log_sig_DF_FC <- function(DF){
  require(ggplot2)
  DFtitle <- deparse(substitute(DF))
  DF <- DF[order(DF$value_1),]
  DF$rank <- (1:length(DF$value_1))
  DF$delta_FPKM <- log(1+abs(DF$value_2 - DF$value_1),2) * (DF$value_2 - DF$value_1)/abs(DF$value_2 - DF$value_1)
  sig_index <- which((DF$significant == "yes" & abs(DF$log) >= 1 & DF$value_1 >= 1)|(DF$significant == "yes" & abs(DF$log) >= 1 & DF$value_2 >= 1)) #same filtering as was applied for network analysis
  notsig_index <- which((DF$significant == "yes" & abs(DF$log) < 1)|(DF$significant == "yes" & DF$value_1 < 1 & DF$value_2 < 1))
  up_sig <- subset(DF[sig_index,],log2.fold_change. > 0)
  up_sig <- round(up_sig$value_2 - up_sig$value_1,2) #list of dFPKMs
  down_sig <- subset(DF[sig_index,],log2.fold_change. < 0)
  down_sig <- round(down_sig$value_2 - down_sig$value_1,2) #list of dFPKMs
  up_notsig <- subset(DF[notsig_index,],log2.fold_change. > 0)
  up_notsig <- round(up_notsig$value_2 - up_notsig$value_1,2) #list of dFPKMs
  down_notsig <- subset(DF[notsig_index,],log2.fold_change. < 0)
  down_notsig <- round(down_notsig$value_2 - down_notsig$value_1,2) #list of dFPKMs
  total.diff <- paste("total DiffExpr(",length(up_sig)+length(down_sig)+length(up_notsig)+length(down_notsig)," genes ; ",sum(up_sig)+sum(down_sig)+sum(up_notsig)+sum(down_notsig)," dFPKM)",sep="")
  total.up.diff <- paste("total UpExpr(",length(up_sig)+length(up_notsig)," genes ; ",sum(up_sig)+sum(up_notsig)," dFPKM)",sep="")
  total.down.diff <- paste("total DownExpr(",length(down_sig)+length(down_notsig)," genes ; ",sum(down_sig)+sum(down_notsig)," dFPKM)",sep="")
  filtered.diff <- paste("total filtered DiffExpr(",length(up_sig)+length(down_sig)," genes ; ",sum(up_sig)+sum(down_sig)," dFPKM)",sep="")
  filtered.up.diff <- paste("Filtered UpExpr(",length(up_sig)," genes ; ",sum(up_sig)," dFPKM)",sep="")
  filtered.down.diff <- paste("Filtered DownExpr(",length(down_sig)," genes ; ",sum(down_sig)," dFPKM)",sep="")
  label_list <- c(total.diff,total.up.diff,total.down.diff,"",filtered.diff,filtered.up.diff,filtered.down.diff)
  
  g = ggplot() +
    geom_point(aes(DF$rank[sig_index],DF$delta_FPKM[sig_index]),size = 1, alpha=0.2,color = "red") +
    geom_point(aes(DF$rank[notsig_index],DF$delta_FPKM[notsig_index]),size = 1, alpha=0.2,color = "black") +
    geom_smooth(color = "red", aes(DF$rank[sig_index],DF$delta_FPKM[sig_index])) +
    geom_smooth(color = "blue", aes(DF$rank[union(sig_index,notsig_index)],DF$delta_FPKM[union(sig_index,notsig_index)])) + #fixed to have q_value cutoff
    geom_hline(yintercept = 0,color = 'black')+
    geom_vline(xintercept = (min(which(DF$value_1 >= 1))),color = 'blue')+ #marks FPKM = 1
    geom_vline(xintercept = (min(which(DF$value_1 >= 10))),color = 'blue')+ #marks FPKM = 10
    geom_vline(xintercept = (min(which(DF$value_1 >= 100))),color = 'blue')+ #marks FPKM = 100
    geom_vline(xintercept = (min(which(DF$value_1 >= 1000))),color = 'blue')+ #marks FPKM = 1000
    geom_vline(xintercept = (min(which(DF$value_1 >= 10000))),color = 'blue')+ #marks FPKM = 10000
    geom_vline(xintercept = (min(which(DF$value_1 >= 100000))),color = 'blue')+ #marks FPKM = 100000
    labs(x = "Rank order: CTRL expression(FPKM)", y = "log2(1+delta_FPKM)", title = DFtitle)+
    annotate("text",x = 4000, y = 17:11, label = label_list)+
    theme(legend.position = "none")+
    theme_bw()
  return(g)
}

##### exclude non-sig genes #######################################used for printed images!!!!!!!!!!!!!!!!!!!!
plot.ranked_expression_log_sig_DF_FC <- function(DF){
  require(ggplot2)
  DFtitle <- deparse(substitute(DF))
  DF <- subset(DF,significant == "yes")
  DF <- DF[order(DF$value_1),]
  DF$rank <- (1:length(DF$value_1))
  DF$delta_FPKM <- log(1+abs(DF$value_2 - DF$value_1),2) * (DF$value_2 - DF$value_1)/abs(DF$value_2 - DF$value_1)
  sig_index <- which((DF$significant == "yes" & abs(DF$log) >= 1 & DF$value_1 >= 1)|(DF$significant == "yes" & abs(DF$log) >= 1 & DF$value_2 >= 1)) #same filtering as was applied for network analysis
  notsig_index <- which((DF$significant == "yes" & abs(DF$log) < 1)|(DF$significant == "yes" & DF$value_1 < 1 & DF$value_2 < 1))
  up_sig <- subset(DF[sig_index,],log2.fold_change. > 0)
  up_sig <- round(up_sig$value_2 - up_sig$value_1,2) #list of dFPKMs
  down_sig <- subset(DF[sig_index,],log2.fold_change. < 0)
  down_sig <- round(down_sig$value_2 - down_sig$value_1,2) #list of dFPKMs
  up_notsig <- subset(DF[notsig_index,],log2.fold_change. > 0)
  up_notsig <- round(up_notsig$value_2 - up_notsig$value_1,2) #list of dFPKMs
  down_notsig <- subset(DF[notsig_index,],log2.fold_change. < 0)
  down_notsig <- round(down_notsig$value_2 - down_notsig$value_1,2) #list of dFPKMs
  total.diff <- paste("total DiffExpr(",length(up_sig)+length(down_sig)+length(up_notsig)+length(down_notsig)," genes ; ",sum(up_sig)+sum(down_sig)+sum(up_notsig)+sum(down_notsig)," dFPKM)",sep="")
  total.up.diff <- paste("total UpExpr(",length(up_sig)+length(up_notsig)," genes ; ",sum(up_sig)+sum(up_notsig)," dFPKM)",sep="")
  total.down.diff <- paste("total DownExpr(",length(down_sig)+length(down_notsig)," genes ; ",sum(down_sig)+sum(down_notsig)," dFPKM)",sep="")
  filtered.diff <- paste("total filtered DiffExpr(",length(up_sig)+length(down_sig)," genes ; ",sum(up_sig)+sum(down_sig)," dFPKM)",sep="")
  filtered.up.diff <- paste("Filtered UpExpr(",length(up_sig)," genes ; ",sum(up_sig)," dFPKM)",sep="")
  filtered.down.diff <- paste("Filtered DownExpr(",length(down_sig)," genes ; ",sum(down_sig)," dFPKM)",sep="")
  label_list <- c(total.diff,total.up.diff,total.down.diff,"",filtered.diff,filtered.up.diff,filtered.down.diff)
  
  g = ggplot() +
    geom_point(aes(DF$rank[sig_index],DF$delta_FPKM[sig_index]),size = 1, alpha=0.2,color = "red") +
    geom_point(aes(DF$rank[notsig_index],DF$delta_FPKM[notsig_index]),size = 1, alpha=0.2,color = "black") +
    geom_smooth(color = "red", aes(DF$rank[sig_index],DF$delta_FPKM[sig_index])) +
    geom_smooth(color = "blue", aes(DF$rank[union(sig_index,notsig_index)],DF$delta_FPKM[union(sig_index,notsig_index)])) + #fixed to have q_value cutoff
    geom_hline(yintercept = 0,color = 'black')+
    geom_vline(xintercept = (min(which(DF$value_1 >= 1))),color = 'blue')+ #marks FPKM = 1
    geom_vline(xintercept = (min(which(DF$value_1 >= 10))),color = 'blue')+ #marks FPKM = 10
    geom_vline(xintercept = (min(which(DF$value_1 >= 100))),color = 'blue')+ #marks FPKM = 100
    geom_vline(xintercept = (min(which(DF$value_1 >= 1000))),color = 'blue')+ #marks FPKM = 1000
    geom_vline(xintercept = (min(which(DF$value_1 >= 10000))),color = 'blue')+ #marks FPKM = 10000
    geom_vline(xintercept = (min(which(DF$value_1 >= 100000))),color = 'blue')+ #marks FPKM = 100000
    labs(x = "Rank order: CTRL expression(FPKM)", y = "log2(1+delta_FPKM)", title = DFtitle)+
    annotate("text",x = 0, y = -7:-13, label = label_list, hjust = 0)+
    theme(legend.position = "none")+
    theme_bw()
  return(g)
}


############### no filtering add dFPKM/gene stat
plot.ranked_expression_log_only.q <- function(DF){
  require(ggplot2)
  DFtitle <- deparse(substitute(DF))
  DF <- subset(DF,significant == "yes") #filter: only q<=0.05 kept
  DF <- DF[order(DF$value_1),]
  DF$rank <- (1:length(DF$value_1))
  DF$delta_FPKM <- log(1+abs(DF$value_2 - DF$value_1),2) * (DF$value_2 - DF$value_1)/abs(DF$value_2 - DF$value_1)
  up <- subset(DF,log2.fold_change. > 0)
  up <- round(up$value_2 - up$value_1,2) #list of dFPKMs
  down <- subset(DF,log2.fold_change. < 0)
  down <- round(down$value_2 - down$value_1,2) #list of dFPKMs
  total.diff <- paste("total DiffExpr(",length(up)+length(down)," genes ; ",sum(up)+sum(down)," dFPKM) ; ",round((sum(up)+sum(down))/(length(up)+length(down)),2)," dFPKM/gene)",sep="")
  total.up.diff <- paste("total UpExpr(",length(up)," genes ; ",sum(up)," dFPKM) ; ",round((sum(up))/(length(up)),2)," dFPKM/gene)",sep="")
  total.down.diff <- paste("total DownExpr(",length(down)," genes ; ",sum(down)," dFPKM) ; ",round((sum(down))/(length(down)),2)," dFPKM/gene)",sep="")
  label_list <- c(total.diff,total.up.diff,total.down.diff)
  
  g = ggplot() +
    geom_point(aes(DF$rank,DF$delta_FPKM),size = 1, alpha=0.2,color = "black") +
    geom_smooth(color = "blue", aes(DF$rank,DF$delta_FPKM)) +
    geom_hline(yintercept = 0,color = 'black')+
    geom_vline(xintercept = (min(which(DF$value_1 >= 1))),color = 'blue')+ #marks FPKM = 1
    geom_vline(xintercept = (min(which(DF$value_1 >= 10))),color = 'blue')+ #marks FPKM = 10
    geom_vline(xintercept = (min(which(DF$value_1 >= 100))),color = 'blue')+ #marks FPKM = 100
    geom_vline(xintercept = (min(which(DF$value_1 >= 1000))),color = 'blue')+ #marks FPKM = 1000
    geom_vline(xintercept = (min(which(DF$value_1 >= 10000))),color = 'blue')+ #marks FPKM = 10000
    geom_vline(xintercept = (min(which(DF$value_1 >= 100000))),color = 'blue')+ #marks FPKM = 100000
    labs(x = "Rank order: CTRL expression(FPKM)", y = "log2(1+delta_FPKM)", title = DFtitle)+
    annotate("text",x = 0, y = -11:-13, label = label_list, hjust = 0)+
    theme(legend.position = "none")+
    theme_bw()
  return(g)
}


#####################################################
#####################################################





#create function that takes diff_expression data ranks by CTRL expression value and plots change in FPKM
#not on log scale, not considering significance
plot.ranked_expression <- function(DFtitle){
  require(ggplot2)
  DF <- get(DFtitle)
  DF <- DF[which(abs(DF$value_2-DF$value_1)>1),]
  DF <- DF[order(DF$value_1),]
  rank <- (1:length(DF$value_1))
  delta_FPKM <- (DF$value_2 - DF$value_1)
  g = ggplot() +
    geom_point(aes(x=rank, y=delta_FPKM), size = 1, alpha=0.3) +
    geom_smooth(aes(x=rank, y=delta_FPKM)) +
    ylim(-1500, 1500) +
    labs(x = "Rank order: CTRL expression(FPKM)", y = "delta_FPKM", title = DFtitle)+
    theme(legend.position = "none")+
    theme_bw()
  return(g)
}

#create function that plots same as above but on log1p scale for delta_FPKM
#all data regardless of significance
plot.ranked_expression_log <- function(DFtitle){
  require(ggplot2)
  DF <- get(DFtitle)
  DF <- DF[order(DF$value_1),]
  DF$rank <- (1:length(DF$value_1))
  DF$delta_FPKM <- log1p(abs(DF$value_2 - DF$value_1)) * (DF$value_2 - DF$value_1)/abs(DF$value_2 - DF$value_1)
  g = ggplot(DF,aes(rank,delta_FPKM)) +
    geom_point(size = 1, alpha=0.2) +
    geom_smooth() +
    geom_hline(yintercept = 0,color = 'red')+
    geom_vline(xintercept = (min(which(DF$value_1 >= 1))),color = 'green')+
    labs(x = "Rank order: CTRL expression(FPKM)", y = "log1p(delta_FPKM)", title = DFtitle)+
    theme(legend.position = "none")+
    theme_bw()
  return(g)
}

#create function that plots same as above but on log1p scale for delta_FPKM
#all data regardless of significance
plot.ranked_expression_log_sig <- function(DFtitle){
  require(ggplot2)
  DF <- get(DFtitle)
  DF <- DF[order(DF$value_1),]
  DF$rank <- (1:length(DF$value_1))
  DF$delta_FPKM <- log1p(abs(DF$value_2 - DF$value_1)) * (DF$value_2 - DF$value_1)/abs(DF$value_2 - DF$value_1)
  sig_index <- which(DF$significant == "yes")
  notsig_index <- which(DF$significant == "no")
  g = ggplot() +
    geom_point(aes(DF$rank[sig_index],DF$delta_FPKM[sig_index]),size = 1, alpha=0.2,color = "red") + #significant genes
    geom_point(aes(DF$rank[notsig_index],DF$delta_FPKM[notsig_index]),size = 1, alpha=0.2,color = "black") + #insignificant genes
    geom_smooth(color = "red", aes(DF$rank[sig_index],DF$delta_FPKM[sig_index])) + # sig data only trendline
    geom_smooth(color = "blue", aes(DF$rank,DF$delta_FPKM)) + #all data trendline
    geom_hline(yintercept = 0,color = 'red')+ #x axis line overlayed on top of data
    geom_vline(xintercept = (min(which(DF$value_1 >= 1))),color = 'green')+ #marks FPKM = 1
    geom_vline(xintercept = (min(which(DF$value_1 >= 10))),color = 'green')+ #marks FPKM = 10
    geom_vline(xintercept = (min(which(DF$value_1 >= 100))),color = 'green')+ #marks FPKM = 100
    geom_vline(xintercept = (min(which(DF$value_1 >= 1000))),color = 'green')+ #marks FPKM = 1000
    geom_vline(xintercept = (min(which(DF$value_1 >= 10000))),color = 'green')+ #marks FPKM = 10000
    geom_vline(xintercept = (min(which(DF$value_1 >= 100000))),color = 'green')+ #marks FPKM = 100000
    labs(x = "Rank order: CTRL expression(FPKM)", y = "log2(1+delta_FPKM)", title = DFtitle)+
    theme(legend.position = "none")+
    theme_bw()
  return(g)
}

#create function that plots on log1p scale for delta_FPKM against actual fpkm values
#all data regardless of significance
plot.actual_expression_log <- function(DFtitle){
  require(ggplot2)
  DF <- get(DFtitle)
  DF$delta_FPKM <- log(1+abs(DF$value_2 - DF$value_1),2) * (DF$value_2 - DF$value_1)/abs(DF$value_2 - DF$value_1)
  g = ggplot(DF,aes(value_1,delta_FPKM)) +
    geom_point(size = 1, alpha=0.2) +
    geom_smooth() +
    geom_hline(yintercept = 0,color = 'red')+
    geom_vline(xintercept = 1,color = 'green')+
    labs(x = "CTRL expression(FPKM)", y = "log2(1+delta_FPKM)", title = DFtitle)+
    scale_x_log10()+
    theme(legend.position = "none")+
    theme_bw()
  return(g)
}

#create function that plots on log1p scale for delta_FPKM against actual fpkm values
#with distinction between significant an not significant
plot.actual_expression_log_sig <- function(DFtitle){
  require(ggplot2)
  DF <- get(DFtitle)
  DF$delta_FPKM <- log(1+abs(DF$value_2 - DF$value_1),2) * (DF$value_2 - DF$value_1)/abs(DF$value_2 - DF$value_1)
  sig_index <- which(DF$significant == "yes")
  notsig_index <- which(DF$significant == "no")
  g = ggplot() +
    geom_point(aes(DF$value_1[sig_index],DF$delta_FPKM[sig_index]),size = 1, alpha=0.2,color = "red") +
    geom_point(aes(DF$value_1[notsig_index],DF$delta_FPKM[notsig_index]),size = 1, alpha=0.2,color = "black") +
    geom_smooth(color = "red", aes(DF$value_1[sig_index],DF$delta_FPKM[sig_index])) +
    geom_smooth(color = "blue", aes(DF$value_1,DF$delta_FPKM)) +
    geom_hline(yintercept = 0,color = 'red')+
    geom_vline(xintercept = 1,color = 'green')+
    labs(x = "CTRL expression(FPKM)", y = "log2(1+delta_FPKM)", title = DFtitle)+
    scale_x_log10()+
    theme(legend.position = "none")+
    theme_bw()
  return(g)
}


##########################################################################################################
##################code to get list of most consistantly downregulated @ G b/t lines###################
##########################################################################################################

##### exclude non-sig genes #######################################used for printed images!!!!!!!!!!!!!!!!!!!!
plot.ranked_expression_log_sig_DF_FC_return.data <- function(DF){
  require(ggplot2)
  DFtitle <- deparse(substitute(DF))
  DF <- subset(DF,significant == "yes")
  DF <- DF[order(DF$value_1),]
  DF$rank <- (length(DF$value_1):1) #in this iteration lowest rank is highest value_1
  DF$delta_FPKM <- (DF$value_2 - DF$value_1)
  sig_index <- which((DF$significant == "yes" & abs(DF$log) >= 1 & DF$value_1 >= 1)|(DF$significant == "yes" & abs(DF$log) >= 1 & DF$value_2 >= 1)) #same filtering as was applied for network analysis
  notsig_index <- which((DF$significant == "yes" & abs(DF$log) < 1)|(DF$significant == "yes" & DF$value_1 < 1 & DF$value_2 < 1))
  return(DF)
}

head(plot.ranked_expression_log_sig_DF_FC_return.data(MA5vMG5))

master.up.down.DF <- function(names.list){ #list of File Names to creat combined and sorted DF
  DF.out <- data.frame(row.names = get(names.list[1])[,"gene"])
  DF.out$sig.freq <-0
  DF.out$average.dFPKM <- 0
  DF.out$average.rank <- 0 
  for(x in names.list){
    DF.in <- get(x)
    DF.in <- plot.ranked_expression_log_sig_DF_FC_return.data(DF.in)
    DF.in$delta <- DF.in$delta*sign(DF.in$log)
    DF.in <- DF.in[DF.in$sig == "yes",c(3,15,16,14)]
    row.names(DF.in) <- DF.in[,1]
    DF.out[row.names(DF.in),"average.rank"] <- DF.out[row.names(DF.in),"average.rank"] + DF.in$rank
    DF.out[row.names(DF.in),"sig.freq"] <- DF.out[row.names(DF.in),"sig.freq"] + 1
    DF.out[row.names(DF.in),"average.dFPKM"] <- DF.out[row.names(DF.in),"average.dFPKM"] + DF.in$delta
    new.col.names <- paste(x,colnames(DF.in),sep = " ")
    DF.out[,new.col.names] <- 0
    DF.out[row.names(DF.in),new.col.names]<-DF.in
  }
  DF.out$average.rank <- DF.out$average.rank/length(names.list)
  DF.out$average.dFPKM <- DF.out$average.dFPKM/length(names.list)
  DF.out <- DF.out[order(-DF.out$sig.freq,-(abs(DF.out$average.dFPKM))),]
  return(DF.out)
}
b <- master.up.down.DF(c("MA5vMG5","FA5vFG5"))
a <- master.up.down.DF(c("MA5vME5","MA5vMG5","FA5vFE5","FA5vFG5"))
