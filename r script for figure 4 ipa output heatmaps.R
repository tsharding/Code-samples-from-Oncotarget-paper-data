# for figure 4 network analysis heat maps. 


#use output from IPA output analysis script copied and pasted (changed to add z-score to GO)


#######################################################################################################################
###################### pasted from other script doc ########changed to add z-score to GO #############
######################################################################################################################

#read in tab delim files
#pathway
#set wd to pathway analysis files
setwd("~/Data after 11_11_15/big exp/RNA-seq/IPA output tab files/pathway analysis")
for(x in list.files()){
  xx <- unlist(strsplit(x,split = ".txt"))
  y <- read.table(x,header = TRUE,sep = "\t",stringsAsFactors = FALSE,skip = 2)
  #get molecules as list instead of comma delim
  y$list <- sapply(y$Molecules, function(x){unlist(strsplit(x,split = ","))})
  y <- y[,c("Ingenuity.Canonical.Pathways","X.log.p.value.","z.score","list")]
  assign(paste(xx,"IPA_pathway",sep = "_"),y)
}
rm(x,xx,y)

#upstream
#set wd to upstream analysis files
setwd("~/Data after 11_11_15/big exp/RNA-seq/IPA output tab files/upstream analysis")
for(x in list.files()){
  xx <- unlist(strsplit(x,split = ".txt"))
  y <- read.table(x,header = TRUE,sep = "\t",stringsAsFactors = FALSE,skip = 2)
  #get molecules as list instead of comma delim
  y$list <- sapply(y$Target.molecules.in.dataset, function(x){unlist(strsplit(x,split = ","))})
  y <- y[,c("Upstream.Regulator","p.value.of.overlap","Activation.z.score","list")]
  assign(paste(xx,"IPA_upstream",sep = "_"),y)
}
rm(x,xx,y)

#GO
#set wd to GO analysis files
setwd("~/Data after 11_11_15/big exp/RNA-seq/IPA output tab files/GO analysis")
for(x in list.files()){
  xx <- unlist(strsplit(x,split = ".txt"))
  y <- read.table(x,header = TRUE,sep = "\t",stringsAsFactors = FALSE,skip = 2,quote = "")
  #get molecules as list instead of comma delim
  y$list <- sapply(y$Molecules, function(x){unlist(strsplit(x,split = ","))})
  print(colnames(y))
  ###################################### changed from origional #######################################################################
  if(length(colnames(y))<9){y[,"Activation.z.score"] <- 0} #if z-score column omitted (all zero NA when exported) add as 0 populated col
  y <- y[,c("Diseases.or.Functions.Annotation","p.Value","Activation.z.score","X","list")]
  ######################################################################################################################################
  assign(paste(xx,"IPA_GO",sep = "_"),y)
}
rm(x,xx,y)

#make lists of df names:
IPA.all <- unlist(strsplit(list.files(),split = ".txt"))
IPA.pathway <- sapply(IPA.all,function(x) {paste(x,"IPA_pathway",sep = "_")})
IPA.upstream <- sapply(IPA.all,function(x) {paste(x,"IPA_upstream",sep = "_")})
IPA.GO <- sapply(IPA.all,function(x) {paste(x,"IPA_GO",sep = "_")})





#######################################################################################################################
################################### generate heatmap ##################################################################
#######################################################################################################################

#pheatmap
install.packages("pheatmap") #installs rcolorbrewer
library("pheatmap")

#accepts list of regulators or pathways in order - plots heat map
#use two color scales (trick into displaying both by hybrid color scale - add 100 to z-scores)
#3 types accepted: "IPA_pathway" ; "IPA_upstream" ; "IPA_GO"
# p-val is column 2 for all three; z-score is column 3 for all three
IPA.heatmap2 <- function(IPA.id.list,type, width = 15, height = 15){ 
  
  DF<-data.frame(row.names = IPA.id.list)#set empty DF to be filled with heatmap data
  names<- c("FA4vFE4","FA5vFE5","FA5vFB5","FA5vFG5","MA4vME4","MA5vME5","MA5vMB5","MA5vMG5")
  for(x in names){
    DF.IPA <- get(paste(x,type,sep = "_"))
    DF.IPA[is.na(DF.IPA)] <- 0 #convert NaN to 0
    if(type == "IPA_upstream" | type == "IPA_GO"){ #subsets for p<0.05
      DF.IPA <- DF.IPA[(DF.IPA[,2] < 0.05),] #exlude non-sig (exclude p>0.05)
    }else{
      DF.IPA <- DF.IPA[(DF.IPA[,2] > (-log10(0.05))),] #exlude non-sig (exclude p>0.05)
    }
    
    row.names(DF.IPA)<- DF.IPA[,1]
    
    if(type == "IPA_upstream" | type == "IPA_GO"){ #converts to -log p_val for upstream and GO to match pathway
      DF[,paste(x,"-log_p_value",sep = " ")] <- -log10(DF.IPA[IPA.id.list,2])
    }else{
      DF[,paste(x,"-log_p_value",sep = " ")] <- DF.IPA[IPA.id.list,2]
    }
    
    DF[is.na(DF)]<-0 # all NA's to 0'
    DF[,paste(x,"activation_z_score",sep = " ")]<-DF.IPA[IPA.id.list,3] + 100
    DF[is.na(DF)]<-0 # all NA's to 0' (not 100!)
  }
  
  DF.ipa.test <<- DF #assign as global variable for troubleshooting
  
  #set color parameters
  hm.breaks <- c(0,seq(-log10(0.05),max(DF[DF<75]),length = 10)) #sets breaks for color scale, total # of breaks = 31 
  #if no negative z scores in entire data frame exclude those breaks/colors
  if(length(DF[DF>75 & DF<100]) != 0){ 
    hm.breaks <- c(hm.breaks,seq(min(DF[DF>75])-0.00001,max(DF[DF<98]),length = 10))
  }else{
    hm.breaks <- c(hm.breaks,99.9999)
  }
  hm.breaks <- c(hm.breaks,seq(min(DF[DF>102]),max(DF[DF>102]),length = 10))
  print("color scale breaks:")
  print(hm.breaks)
  
  col1 <- (colorRampPalette(c("black",'red'))(9+4))[5:13]
  col2 <- (colorRampPalette(c("purple",'grey'))(9+4))[1:9]
  col3 <- (colorRampPalette(c("grey",'orange'))(9+4))[5:13]
  hm.colors <- c("white",col1,"blue") #should be 1 less than hm.breaks; white = no enrichment
  #if no negative z scores in entire data frame exclude those breaks/colors
  if(length(DF[DF>75 & DF<100]) != 0){
    hm.colors <- c(hm.colors,col2)
  }
  hm.colors <- c(hm.colors,"white",col3)
  
  #make annotation data frame
  hm.col.annotation <- data.frame(Cell.Line <- factor(rep(c("FLAM76","MMM1"),each = 8)),
                                  Day <- factor(rep(c("Day 4","Day 4",rep("Day 5.5",6)),2)),
                                  EPZ <- factor(rep(c("+","+","+","+","-","-","+","+"),2)),
                                  PAN <- factor(rep(c("-","-","-","-","+","+","+","+"),2)),
                                  row.names = colnames(DF))
  colnames(hm.col.annotation)<-c("Cell Line","Sampling Day     ","EPZ6438","Panobinostat")
  #specify colors for annotation
  hm.col.annotation.colors <- list("Cell Line" = c(FLAM76 = "gold",MMM1 ="darkorange2"),
                                   "Sampling Day     " = c("Day 4" = "royalblue","Day 5.5" = "navy"),
                                   "EPZ6438" = c("+" = "black","-" = "white"),
                                   "Panobinostat" = c("+" = "black","-" = "white"))
  
  #generate heatmap
  if(height < 1) {hm.ylab <- FALSE}else{hm.ylab <- TRUE} #dont display genes if cells are too short
  pheatmap(DF,
           color = hm.colors,
           breaks = hm.breaks,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           gaps_col = c(2,2,2,4,6,8,8,8,8,8,8,10,10,10,12,14),#by day
           cellwidth = width,
           cellheight = height,
           treeheight_row = 0,
           annotation_col = hm.col.annotation,
           annotation_colors = hm.col.annotation.colors,
           show_colnames = FALSE,
           show_rownames = hm.ylab,
           legend_labels = c("log2 Fold change"),
           fontsize_row = height
  )
}

query <- unique(c(MA5vMG5_IPA_GO$Diseases.or.Functions.Annotation[1:10],
                  FA5vFG5_IPA_GO$Diseases.or.Functions.Annotation[1:10],
                  MA5vMB5_IPA_GO$Diseases.or.Functions.Annotation[1:10],
                  FA5vFB5_IPA_GO$Diseases.or.Functions.Annotation[1:10],
                  MA5vME5_IPA_GO$Diseases.or.Functions.Annotation[1:10],
                  FA5vFE5_IPA_GO$Diseases.or.Functions.Annotation[1:10]))

IPA.heatmap2(query,type = "IPA_GO",height = 10,width = 10)

query <- unique(c(MA5vMG5_IPA_pathway$Ingenuity.Canonical.Pathways[1:10],
                  FA5vFG5_IPA_pathway$Ingenuity.Canonical.Pathways[1:10],
                  MA5vMB5_IPA_pathway$Ingenuity.Canonical.Pathways[1:10],
                  FA5vFB5_IPA_pathway$Ingenuity.Canonical.Pathways[1:10],
                  MA5vME5_IPA_pathway$Ingenuity.Canonical.Pathways[1:10],
                  FA5vFE5_IPA_pathway$Ingenuity.Canonical.Pathways[1:10]))

IPA.heatmap2(query,type = "IPA_pathway",height = 10,width = 10)

query <- unique(c(MA5vMG5_IPA_upstream$Upstream.Regulator[1:20],
                  FA5vFG5_IPA_upstream$Upstream.Regulator[1:20],
                  MA5vMB5_IPA_upstream$Upstream.Regulator[1:20],
                  FA5vFB5_IPA_upstream$Upstream.Regulator[1:20],
                  MA5vME5_IPA_upstream$Upstream.Regulator[1:20],
                  FA5vFE5_IPA_upstream$Upstream.Regulator[1:20]))

IPA.heatmap2(query,type = "IPA_upstream",height = 8,width = 8)
IPA.heatmap.legend(1,query,type = "IPA_upstream",height = 18,width = 8)
IPA.heatmap.legend(2,query,type = "IPA_upstream",height = 18,width = 8)


####################################################################################################
###################################dummy heatmaps for legends#######################################
###################################################################################################
#same as above but with different breaks and colors (z-score scale normalized around zero)

IPA.heatmap.legend <- function(legend_number,IPA.id.list,type, width = 15, height = 15){ 
  #accepts 1 or 2 for legend number to return legend for p_val or z_score respectively
  DF<-data.frame(row.names = IPA.id.list)#set empty DF to be filled with heatmap data
  names<- c("FA4vFE4","FA5vFE5","FA5vFB5","FA5vFG5","MA4vME4","MA5vME5","MA5vMB5","MA5vMG5")
  for(x in names){
    DF.IPA <- get(paste(x,type,sep = "_"))
    DF.IPA[is.na(DF.IPA)] <- 0 #convert NaN to 0
    if(type == "IPA_upstream" | type == "IPA_GO"){ #subsets for p<0.05
      DF.IPA <- DF.IPA[(DF.IPA[,2] < 0.05),] #exlude non-sig (exclude p>0.05)
    }else{
      DF.IPA <- DF.IPA[(DF.IPA[,2] > (-log10(0.05))),] #exlude non-sig (exclude p>0.05)
    }
    
    row.names(DF.IPA)<- DF.IPA[,1]
    
    if(type == "IPA_upstream" | type == "IPA_GO"){ #converts to -log p_val for upstream and GO to match pathway
      DF[,paste(x,"-log_p_value",sep = " ")] <- -log10(DF.IPA[IPA.id.list,2])
    }else{
      DF[,paste(x,"-log_p_value",sep = " ")] <- DF.IPA[IPA.id.list,2]
    }
    
    DF[is.na(DF)]<-0 # all NA's to 0'
    DF[,paste(x,"activation_z_score",sep = " ")]<-DF.IPA[IPA.id.list,3] + 100
    DF[is.na(DF)]<-0 # all NA's to 0' (not 100!)
  }
  
  
  #set breaks and colors depending on legend_nuber arg  
  if(legend_number == 1){
    hm.breaks <- c(0,seq(-log10(0.05),max(DF[DF<75]),length = 10))
    col1 <- (colorRampPalette(c("black",'red'))(9+4))[5:13]
    hm.colors <- c("white",col1)
  } 
  
  if(legend_number == 2){
    
    if(length(DF[DF>75 & DF<100]) != 0){ 
      hm.breaks <- (seq(min(DF[DF>75])-0.00001,max(DF[DF<98]),length = 10))
    }else{
      hm.breaks <- c(hm.breaks,99.9999)
    }
    hm.breaks <- c(hm.breaks,seq(min(DF[DF>102]),max(DF[DF>102]),length = 10)) - 100 #corrects back to zero centered scale
    
    col2 <- (colorRampPalette(c("purple",'grey'))(9+4))[1:9]
    col3 <- (colorRampPalette(c("grey",'orange'))(9+4))[5:13]
    if(length(DF[DF>75 & DF<100]) != 0){
      hm.colors <- (col2)
    }
    hm.colors <- c(hm.colors,"white",col3)
  } 
  
  
  #make annotation data frame
  hm.col.annotation <- data.frame(Cell.Line <- factor(rep(c("FLAM76","MMM1"),each = 8)),
                                  Day <- factor(rep(c("Day 4","Day 4",rep("Day 5.5",6)),2)),
                                  EPZ <- factor(rep(c("+","+","+","+","-","-","+","+"),2)),
                                  PAN <- factor(rep(c("-","-","-","-","+","+","+","+"),2)),
                                  row.names = colnames(DF))
  colnames(hm.col.annotation)<-c("Cell Line","Sampling Day     ","EPZ6438","Panobinostat")
  #specify colors for annotation
  hm.col.annotation.colors <- list("Cell Line" = c(FLAM76 = "gold",MMM1 ="darkorange2"),
                                   "Sampling Day     " = c("Day 4" = "royalblue","Day 5.5" = "navy"),
                                   "EPZ6438" = c("+" = "black","-" = "white"),
                                   "Panobinostat" = c("+" = "black","-" = "white"))
  
  #generate heatmap
  if(height < 1) {hm.ylab <- FALSE}else{hm.ylab <- TRUE} #dont display genes if cells are too short
  pheatmap(DF[1:10,],
           color = hm.colors,
           breaks = hm.breaks,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           gaps_col = c(2,2,2,4,6,8,8,8,8,8,8,10,10,10,12,14),#by day
           cellwidth = width,
           cellheight = height,
           treeheight_row = 0,
           show_colnames = FALSE,
           show_rownames = FALSE,
           fontsize_row = height
  )
}

####################################################################################################
################################### function to 'cluster' highly enriched results #################
###################################################################################################

#used below: list of variable names generated in IPA output analysis script
  #IPA.pathway 
  #IPA.upstream 
  #IPA.GO

#get union of all results for a certain type 
#pathway
results.union <- list()
for(x in IPA.pathway){
  results.union <- unlist(c(results.union,(get(x))[,"Ingenuity.Canonical.Pathways"]))
}
IPA.pathway.results.union <- unique(unlist(results.union))
#upstream
results.union <- list()
for(x in IPA.upstream){
  results.union <- unlist(c(results.union,(get(x))[,"Upstream.Regulator"]))
}
IPA.upstream.results.union <-unique(unlist(results.union))
#GO
  results.union <- list()
for(x in IPA.GO){
  results.union <- unlist(c(results.union,(get(x))[,"Diseases.or.Functions.Annotation"]))
}
IPA.GO.results.union <- unique(unlist(results.union))

rm(x,results.union)  
  
IPA.cluster <- function(type,line){
  #returns list of ipa results to input into heatmap
  #type -> 3 types accepted: "pathway" ; "upstream" ; "GO"
  #line: 1 = FLAM76, 2 = MMM1, 3 = both - scope of enrichment, 4 = AvG
  results.list <- get(paste("IPA",type,"results.union", sep = "."))
  DF<-data.frame(row.names = results.list)#set empty DF to be filled with heatmap data
  names<- c("FA4vFE4","FA5vFE5","FA5vFB5","FA5vFG5","MA4vME4","MA5vME5","MA5vMB5","MA5vMG5")
  type2<- paste("IPA",type,sep = "_")
  for(x in names){
    DF.IPA <- get(paste(x,type2,sep = "_"))
    DF.IPA[is.na(DF.IPA)] <- 0 #convert NaN to 0
    if(type2 == "IPA_upstream" | type2 == "IPA_GO"){ #subsets for p<0.05
      DF.IPA <- DF.IPA[(DF.IPA[,2] < 0.05),] #exlude non-sig (exclude p>0.05)
    }else{
      DF.IPA <- DF.IPA[(DF.IPA[,2] > (-log10(0.05))),] #exlude non-sig (exclude p>0.05)
    }
    
    row.names(DF.IPA)<- DF.IPA[,1]
    
    if(type2 == "IPA_upstream" | type2 == "IPA_GO"){ #converts to -log p_val for upstream and GO to match pathway
      DF[,paste(x,"-log_p_value",sep = " ")] <- -log10(DF.IPA[results.list,2])
    }else{
      DF[,paste(x,"-log_p_value",sep = " ")] <- DF.IPA[results.list,2]
    }
    
    DF[is.na(DF)]<-0 # all NA's to 0'
    DF[,paste(x,"activation_z_score",sep = " ")]<-DF.IPA[results.list,3] #not scaled around 100 as in heatmap
    DF[is.na(DF)]<-0 # all NA's to 0' (not 100!)
  }
  DF$p_value_F_total <- (DF[,1] != 0)+(DF[,3] != 0)+(DF[,5] != 0)+(DF[,7] != 0)
  DF$p_value_M_total <- (DF[,9] != 0)+(DF[,11] != 0)+(DF[,13] != 0)+(DF[,15] != 0)
  DF$p_value_total <- DF$p_value_F_total + DF$p_value_M_total
  DF$p_value_F_average <- (DF[,1]+DF[,3]+DF[,5]+DF[,7])/4
  DF$p_value_M_average <- (DF[,9]+DF[,11]+DF[,13]+DF[,15])/4
  DF$p_value_average <- (DF[,1]+DF[,3]+DF[,5]+DF[,7]+DF[,9]+DF[,11]+DF[,13]+DF[,15])/8
  DF$z_score_F_total <-  (DF[,2] != 0)+(DF[,4] != 0)+(DF[,6] != 0)+(DF[,8] != 0)
  DF$z_score_M_total <-  (DF[,10] != 0)+(DF[,12] != 0)+(DF[,14] != 0)+(DF[,16] != 0)
  DF$z_score_total <- DF$z_score_F_total + DF$z_score_M_total
  DF$G_only <- (DF[,7] != 0)+(DF[,15] != 0)
  #order by FLAM 
  if(line == 1){
    DF <- DF[order(-DF$p_value_F_total,DF$p_value_M_total),]
  }
  #order by MMM1
  if(line == 2){
    DF <- DF[order(-DF$p_value_M_total,DF$p_value_F_total),]
  }
  #order by both
  if(line == 3){
    DF <- DF[order(-DF$p_value_total,-DF$p_value_average),]
  }
  #order for G only
  if(line == 4){
    DF <- DF[order(-DF$G_only,DF$p_value_total),]
  }
  
  DF.ipa.cluster.test <<- DF
  return(row.names(DF))
}

head(IPA.cluster("GO",3))

############################################################
IPA.heatmap2(IPA.cluster("upstream",3)[1:100],type = "IPA_upstream",height = 6,width = 8)

IPA.final.upstream.list <-c(
  "EZH2",
  "Hdac",
  "HDAC1",
  "HDAC2",
  "CCND1",
  "TP53",
  "MYC",
  "BTK",
  "MAPK1",
  "ERBB2",
  "TNF",
  "ATM",
  "Jnk",
  "BRD4",
  "P38 MAPK",
  "EIF2AK2",
  "SMARCA4",
  "SMARCB1",
  "Interferon alpha",
  "IFNG",
  "STAT1",
  "IRF4",
  "SOCS1",
  "SOCS3",
  "estrogen receptor",
  "ESR1",
  "BCL6",
  "PTEN",
  "AURKB",
  "HIF1A",
  "TGFB1",
  "NFkB (complex)"
)
IPA.final.upstream.list[which(!(IPA.final.upstream.list %in% IPA.upstream.results.union))]
IPA.heatmap2(IPA.final.upstream.list,type = "IPA_upstream",height = 15,width = 15)
IPA.heatmap.legend(1,IPA.final.upstream.list,type = "IPA_upstream",height = 18,width = 8)
IPA.heatmap.legend(2,IPA.final.upstream.list,type = "IPA_upstream",height = 18,width = 8)


IPA.final.GO.list <-c(
  "multiple myeloma",
  "cell death",
  "cancer",
  "invasion of tumor cell lines",
  "colony formation of tumor cell lines",
  "cell movement",
  "migration of cells",
  "antimicrobial response",
  "non-Hodgkin disease",
  "hematologic cancer",
  "adhesion of tumor cell lines",
  "colony formation of cells",
  "B-cell leukemia",
  "tumorigenesis of tissue",
  "cytolysis",
  "cell movement of leukocytes",
  "cell death of immune cells",
  "adhesion of immune cells",
  "replication of RNA virus",
  "apoptosis of leukemia cell lines",
  "cell death of leukemia cell lines",
  "degradation of DNA",
  "formation of actin stress fibers",
  "arrest in interphase",
  "arrest in G1 phase of tumor cell lines",
  "cell cycle progression",
  "expression of RNA"
)
IPA.final.GO.list[which(!(IPA.final.GO.list %in% IPA.GO.results.union))]
IPA.heatmap2(IPA.final.GO.list,type = "IPA_GO",height = 15,width = 15)
IPA.heatmap.legend(1,IPA.final.GO.list,type = "IPA_GO",height = 18,width = 8)
IPA.heatmap.legend(2,IPA.final.GO.list,type = "IPA_GO",height = 18,width = 8)

IPA.final.pathway.list <-c(
  "Interferon Signaling",
  "B Cell Development",
  "Antigen Presentation Pathway",
  "Molecular Mechanisms of Cancer",
  "Calcium-induced T Lymphocyte Apoptosis",
  "IL-4 Signaling",
  "Th1 and Th2 Activation Pathway",
  "Integrin Signaling",
  "Phagosome Formation",
  "Axonal Guidance Signaling",
  "Epithelial Adherens Junction Signaling",
  "IL-8 Signaling",
  "Granzyme A Signaling",
  "Estrogen-mediated S-phase Entry",
  "Cyclins and Cell Cycle Regulation",
  "Cell Cycle Control of Chromosomal Replication",
  "Fatty Acid β-oxidation I"
)
IPA.final.pathway.list[which(!(IPA.final.pathway.list %in% IPA.pathway.results.union))]
IPA.heatmap2(IPA.final.pathway.list,type = "IPA_pathway",height = 15,width = 15)
IPA.heatmap.legend(1,IPA.final.pathway.list,type = "IPA_pathway",height = 18,width = 8)
IPA.heatmap.legend(2,IPA.final.pathway.list,type = "IPA_pathway",height = 18,width = 8)

