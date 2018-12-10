# Unsupervised learning methods
# PCA approach
rm(list=ls())
gene <- read.csv("C:/Users/User 1/Desktop/gene.csv")
pheno <- read.csv("C:/Users/User 1/Desktop/pheno.csv")
head(gene)
dim(gene)
head(pheno)
length(pheno)
wang.pheno=pheno$cell.type
gene$Gene <- gene$Gene
gene$Gene <- NULL
head(gene)
dim(gene)
wang.data=t(gene)

dim(wang.data)
head(wang.data)
table(wang.pheno)

pr.out2 =prcomp (wang.data, scale=FALSE)

Cols=function(vec){
  cols=rainbow (length(unique(vec )))
  return (cols[as.numeric(as.factor (vec))])
}

par(mfrow =c(1,2))
plot(pr.out2$x[,1:2], col =Cols(wang.pheno), pch =19, xlab ="Z1",ylab="Z2")
plot(pr.out2$x[,c(1,3)],col =Cols(wang.pheno), pch =19, xlab ="Z1",ylab="Z3")

summary (pr.out2)

plot(pr.out2)
summary(pr.out2)$importance[2,]
summary(pr.out2)$importance[3,]
# Scree plot
pve = 100* pr.out2$sdev ^2/ sum(pr.out2$sdev ^2)
par(mfrow =c(1,2))
plot(pve , type ="o", ylab="PVE ", xlab="Principal Component", col ="blue")
plot(cumsum (pve ), type="o", ylab ="Cumulative PVE", xlab="Principal Component ", col ="brown3 ")

