rm(list = ls())
library(glmnet)
library(plyr)
library(caret)
library(ggplot2)
library(Hmisc)
library(survival)
library(tidyverse)
gene<-read.csv("gene.csv")
gene<-as.data.frame(t(gene))
gene<-gene[c(1:25),3]
gene<-as.data.frame(gene)
gene<-gene$gene
colnames(gene)<-"gene"
expr<-read.csv("teain_test.csv",header = T,row.names = 1)
set.seed(12345)
expr$os_time<-expr$os_time/365
expr3<-na.omit(expr3)
expr3<-filter(expr3,expr3$os_time >= 0)
expr3[]<-lapply(expr3,as.numeric)
expr2<-expr3[,c(3:19)]
expr1<-expr3[,c(1:2)]
expr2<-log2(expr2+1)
expr3<-cbind(expr1,expr2)

write.csv(expr3,file = "dataforlasso.csv")
expr3<-read.csv("dataforlasso.csv",header = T,row.names = 1)
expr<-expr3
x = as.matrix(expr[,gene$gene])
y = data.matrix(Surv(expr$os_time,expr$event))
fit = glmnet(x,y,family = "binomial",nlambda=100,alpha=1)
plot(fit,xvar="lambda",label=T)
rm(weights)
y <- expr$event
cvfit=cv.glmnet(x,y,weights = NULL)
plot(cvfit)
coef<-coef(fit,s=cvfit$lambda.min)
index<-which(coef!=0)
actcoef<-coef[index]
lassoGene<-rownames(coef)[index]
genecoef<- cbind(Gene= lassoGene,Coef=actcoef)
genecoef<-genecoef[-1,]
genecoef<-as.data.frame(genecoef)
write.csv(genecoef,file = "genecoef.csv")
trainFInalGeneExp=expr[,genecoef$Gene]
myFun = function(x){crossprod(as.numeric(x),actcoef)}
actcoef <- actcoef[-1]
Trainscore<-apply(trainFInalGeneExp,1,myFun)
lassoGene<- lassoGene[-1]
outcol<-c("event","os_time",lassoGene)
Risk<-as.vector(ifelse(Trainscore>median(Trainscore),"High","Low"))
outtable<-cbind(expr[,outcol],riskScore=as.vector(Trainscore),Risk)
write.csv(outtable,file = "internal_risK_group.csv")