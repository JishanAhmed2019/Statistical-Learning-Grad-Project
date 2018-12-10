# Unsupervised learning methods
# PCA approach
rm(list=ls())
library(ISLR)
# Case Study III: genomic data

# NCI60 cancer cell line microarray data, 
# 6,830 gene expression measurements on 64 cancer cell lines
#head(NCI60)
gene <- read.csv("C:/Users/User 1/Desktop/gene.csv")
pheno <- read.csv("C:/Users/User 1/Desktop/pheno.csv")
head(gene)
dim(gene)
head(pheno)
pheno
nci.labs=pheno$cell.type
gene$Gene <- gene$Gene
gene$Gene <- NULL
head(gene)
dim(gene)
nci.data=t(gene)

dim(nci.data)
head(nci.data)
table(nci.labs)

pr.out2 =prcomp (nci.data, scale=FALSE)

Cols=function(vec){
  cols=rainbow (length(unique(vec )))
  return (cols[as.numeric(as.factor (vec))])
}

par(mfrow =c(1,2))
plot(pr.out2$x[,1:2], col =Cols(nci.labs), pch =19, xlab ="Z1",ylab="Z2")
plot(pr.out2$x[,c(1,3)],col =Cols(nci.labs), pch =19, xlab ="Z1",ylab="Z3")

#summary (pr.out2)

plot(pr.out2)

# Scree plot
pve = 100* pr.out2$sdev ^2/ sum(pr.out2$sdev ^2)
par(mfrow =c(1,2))
plot(pve , type ="o", ylab="PVE ", xlab="Principal Component", col ="blue")
plot(cumsum (pve ), type="o", ylab ="Cumulative PVE", xlab="Principal Component ", col ="brown3 ")

