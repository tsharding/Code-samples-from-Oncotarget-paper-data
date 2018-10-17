#convert entrez to symbol
#https://www.r-bloggers.com/converting-gene-names-in-r-with-annotationdbi/

source("http://bioconductor.org/biocLite.R")
biocLite("AnnotationDbi")
biocLite("org.Hs.eg.db")
library("AnnotationDbi")#use this package to convert entrez IDs to gen symbols
library("org.Hs.eg.db")


mapIds(org.Hs.eg.db, keys=y[,1], column="SYMBOL", keytype="ENTREZID", multiVals="first")