# venn diagrams for figure 3

install.packages('VennDiagram')
library(VennDiagram)

VenDi3 <- function(DF1,DF2,DF3,FCthresh){ #accepts unfiltered data.frames (3) and FC threshold (not log2!) -> outputs vendiagram
  DF1_filtered <- (filter.df.data(DF1,FCthresh))[,"gene"]
  DF2_filtered <- (filter.df.data(DF2,FCthresh))[,"gene"]
  DF3_filtered <- (filter.df.data(DF3,FCthresh))[,"gene"]
  a1 <- length(DF1_filtered)
  print(a1)
  a2 <- length(DF2_filtered)
  print(a2)
  a3 <- length(DF3_filtered)
  print(a3)
  nn12 <- length(intersect(DF1_filtered,DF2_filtered))
  print(nn12)
  nn23 <- length(intersect(DF2_filtered,DF3_filtered))
  print(nn23)
  nn13 <- length(intersect(DF1_filtered,DF3_filtered))
  print(nn13)
  nn123 <- length(intersect(intersect(DF1_filtered,DF2_filtered),DF3_filtered))
  print(nn123)
  
  
  grid.newpage()
  draw.triple.venn(area1 = a1, area2 = a2, area3 = a3, n12 = nn12, n23 = nn23, n13 = nn13, 
                   n123 = nn123, category = c("EPZ-6438", "Panobinostat", "Combo"),  
                   fill = c("blue", "yellow", "green"),cex = 3,cat.cex = 2.5,cat.pos = c(-10,10,180),cat.dist = c(.1,.1,.1), 
                   euler.d = FALSE,scaled = FALSE,alpha = .4
                   )
}

VenDi3(MA5vME5,MA5vMB5,MA5vMG5,2)
VenDi3(FA5vFE5,FA5vFB5,FA5vFG5,2)

Find # shared b/t unique to combo
UTC_shared <- function(DF1,DF2,DF3,DF4,DF5,DF6,FCthresh){ #accepts unfiltered data.frames (6) and FC threshold (not log2!) -> outputs # of shared unique to combo genes
  DF1_filtered <- (filter.df.data(DF1,FCthresh))[,"gene"]
  print(length(DF1_filtered))
  DF2_filtered <- (filter.df.data(DF2,FCthresh))[,"gene"]
  print(length(DF2_filtered))
  DF3_filtered <- (filter.df.data(DF3,FCthresh))[,"gene"]
  print(length(DF3_filtered))
  DF4_filtered <- (filter.df.data(DF4,FCthresh))[,"gene"]
  print(length(DF4_filtered))
  DF5_filtered <- (filter.df.data(DF5,FCthresh))[,"gene"]
  print(length(DF5_filtered))
  DF6_filtered <- (filter.df.data(DF6,FCthresh))[,"gene"]
  print(length(DF6_filtered))
  #Flam
  UTC_F <- DF3_filtered[!(DF3_filtered %in% union(DF1_filtered,DF2_filtered))]
  print(length(UTC_F))
  #MMM1
  UTC_M <- DF6_filtered[!(DF6_filtered %in% union(DF4_filtered,DF5_filtered))]
  print(length(UTC_M))
  print("shared UTC:")
  print(length(intersect(UTC_F,UTC_M)))
}


UTC_shared(FA5vFE5,FA5vFB5,FA5vFG5,MA5vME5,MA5vMB5,MA5vMG5,2)




