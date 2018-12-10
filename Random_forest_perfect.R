# Clustering methods
rm(list=ls())
# Case Study III: genomic data

# NCI60 cancer cell line microarray data, 
# 6,830 gene expression measurements on 64 cancer cell lines
#head(NCI60)
gene <- read.csv("C:/Users/User 1/Desktop/gene.csv")
pheno <- read.csv("C:/Users/User 1/Desktop/pheno.csv")
head(gene)
dim(gene)
head(pheno)
length(pheno)
wang.pheno=pheno$cell.type
head(wang.pheno)
gene_Gene <- gene[,-1]
head(gene_Gene)
dim(gene_Gene)
wang.data=t(gene_Gene)

dim(wang.data)
colnames(wang.data)
str(wang.data)
table(wang.pheno)

pr.out2 =prcomp (wang.data, scale=FALSE)
names(pr.out2)
pr.out2$rotation
dim(pr.out2$x)
#biplot(pr.out2, scale = 0)
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
dimension.data <- data.frame(wang.pheno, pr.out2$x)
pca.data <- dimension.data[,1:11]
dim(pca.data)
str(pca.data)
#library(rpart)
#rpart.model <- rpart(wang.pheno ~ .,data = train.data, method = "anova")
#rpart.model
#plot(rpart.model)
set.seed(81)
library(MASS)
library(rpart)
library(randomForest)
rf.model <- randomForest(wang.pheno ~ ., data=pca.data, importance=TRUE)
rf.model
summary(rf.model)
plot(rf.model)
colnames(pca.data)
#Prediction on training data
predicted.values <- predict(rf.model, pca.data[1:11])
d_pca <- table(predicted.values, pca.data$wang.pheno)
print(d_pca)
library(caret)
confusionMatrix(predicted.values,pca.data$wang.pheno)

