# Clustering methods
rm(list=ls())
# Case Study III: genomic data

# NCI60 cancer cell line microarray data, 
# 6,830 gene expression measurements on 64 cancer cell lines
#head(NCI60)
gene <- read.csv("C:/Users/User 1/Desktop/gene.csv")
pheno <- read.csv("C:/Users/User 1/Desktop/pheno.csv")
head(gene)
head(pheno)
pheno
wang.pheno=pheno$cell.type
gene$Gene <- gene$Gene
gene$Gene <- NULL
head(gene)
dim(gene)
wang.data=t(gene)
head(wang.data)
# Hierarchical clustering
dim(wang.data)

# Clustering the Observations of the NCI60 Data

#nci.labs=NCI60$labs
#nci.data=NCI60$data

dim(wang.data)
head(wang.data)
dim(wang.data)
# Standardize data to have mean 0 and stdev 1 
sd.data=scale(wang.data)
# Create distance matrix
data.dist=dist(sd.data)

hc.single =hclust(data.dist, method ="single")
hc.complete =hclust(data.dist, method ="complete")
hc.average =hclust(data.dist, method ="average")

par(mfrow =c(3,1))
plot(hc.single, labels =wang.pheno, main="Single Linkage", xlab="",sub ="", ylab="" )
plot(hc.complete ,labels =wang.pheno, main ="Complete Linkage", xlab="", sub ="", ylab="" )
plot(hc.average , labels =wang.pheno, main ="Average Linkage", xlab="", sub ="", ylab="" )

#hc.out
hc.out =hclust(dist(sd.data))
hc.clusters =cutree (hc.out ,3)
table(hc.clusters ,wang.pheno)

par(mfrow =c(1,1))
plot(hc.out , labels =wang.pheno)
abline (h=371, col ="red")

# K means clustering on NCI60 data
any(is.na(sd.data))
sd.data
head(model.matrix(~.+0, sd.data))
set.seed(2)
km.out =kmeans(sd.data, 2, nstart =20)
km.clusters =km.out$cluster

plot(sd.data, col =(km.out$cluster), main="K-Means Clustering
     Results with K=2", xlab ="", ylab="", pch =20, cex =2)

# Compare hierarchical and K-means clustering

table(km.clusters,hc.clusters)

# To get better results run clutering on first few PC

pr.out = prcomp(wang.data, scale=TRUE)

hc.out = hclust(dist(pr.out$x[,1:7]))
plot(hc.out, labels= wang.pheno, main="Hier. Clust. on First seven PC")
table(cutree(hc.out,4),wang.pheno)
