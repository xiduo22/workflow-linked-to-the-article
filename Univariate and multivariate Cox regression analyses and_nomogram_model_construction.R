library(tibble)
library(tidyverse)
a325<-read.table("CGGA.mRNAseq_325.RSEM-genes.20200506.txt",header = T,sep = "\t",row.names = 1)
b693<-read.table("CGGA.mRNAseq_693.RSEM-genes.20200506.txt",header = T,sep = "\t",row.names = 1)
clini_325<-read.table("CGGA.mRNAseq_325_clinical.20200506.txt",header = T,sep = "\t")
clini_693<-read.table("CGGA.mRNAseq_693_clinical.20200506.txt",header = T,sep = "\t")
fix(clini_325)
clini_merge<-rbind(clini_325,clini_693)
gene<-read.csv("genecoef.csv",header = T,row.names = 1)
a325<-a325[gene$Gene,]
b693<-b693[gene$Gene,]
expr<-cbind(a325,b693)
expr<-as.data.frame(t(expr))
expr<-log2(expr+1)
trainFInalGeneExp_all<-expr[,gene$Gene]
myFun = function(x){crossprod(as.numeric(x),gene$Coef)}
Trainscore<-apply(trainFInalGeneExp_all,1,myFun)
outcol<-c(gene$Gene)
Risk<-as.vector(ifelse(Trainscore>median(Trainscore),"High","Low"))
outtableall<-cbind(expr[,outcol],riskScore=as.vector(Trainscore),Risk)
outtableall<-rownames_to_column(outtableall,var = "CGGA_ID")
tb<-inner_join(clini_merge,outtableall,by="CGGA_ID")

tb1<-tb
tb1 <- tb1[tb1$event %in% c(0, 1), ]
tb1 <- subset(tb1, !is.na(PRS_type))
library(survival)
library(forestplot)
library(rms)
write.csv(tb1,file = "cgga .csv")
rm(list = ls())
tb1<-read.csv("cgga .csv",header = T,row.names = 1)
tb1<- na.omit(tb1)
tb1 <- tb1 %>%
mutate(Histology = ifelse(Histology %in% c("GBM", "rGBM"), "GBM", "no_GBM"))
tb1 <- tb1 %>%
mutate(Grade = ifelse(Grade %in% c("WHO III", "WHO IV"), "High", "Low"))
tb1 <- tb1 %>%
  mutate(PRS_type = ifelse(PRS_type %in% c("Primary" ), "Primary", "no_Primary"))

tb1$Gender <- trimws(tb1$Gender)
tb1$Gender <- tolower(tb1$Gender)  

tb1$Gender <- factor(tb1$Gender, levels = c("male", "female"))

tb1$PRS_type <- factor(tb1$PRS_type, levels = c("Primary", "no_Primary"))
tb1$Histology <- factor(tb1$Histology, levels = c( "no_GBM","GBM"))
tb1$Grade <- factor(tb1$Grade, levels = c("Low","High"))

tb1$IDH_mutation_status <- factor(tb1$IDH_mutation_status, levels = c("Wildtype","Mutant"))
tb1$X1p19q_codeletion_status <- factor(tb1$X1p19q_codeletion_status, levels = c("Non-codel", "Codel"))
tb1$Gender <- factor(tb1$Gender, levels = c("male", "female"))
res.cox <- coxph(Surv(OS, event) ~ Age, data = tb1)
res.cox
tb1$riskScore

tb1 <- tb1 %>%
  mutate(Age = ifelse(Age >= 60, "Old", "Young"))

tb1$Age <- factor(tb1$Age, levels = c("Young", "Old"))

covariates <- c("riskScore", "PRS_type", "Histology", "Grade ", "Gender", "Age", "Radio_status", "Chemo_status", "IDH_mutation_status", "X1p19q_codeletion_status")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS, event)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = tb1)})


univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         
                         HR <-signif(x$coef[2], digits=2);
                         
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
write.csv(res,file = "genes_uncox.csv")

tb1$PRS_type <- relevel(tb1$PRS_type, ref = "Primary")
res.cox <- coxph(Surv(OS, event) ~ riskScore+PRS_type+Histology+Grade+Gender+Age+Radio_status+Chemo_status+IDH_mutation_status+X1p19q_codeletion_status, data =  tb1)
res.cox
x <- summary(res.cox)
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)

multi_res=data.frame(p.value=pvalue,
                     HR=paste(HR," (",low,"-",high,")",sep=""),
                     stringsAsFactors = F
)
multi_res
write.csv(multi_res,file = "multicox.csv")

library(VIM)
library(nomogramFormula)
library(rmda)
library(pROC)
library(rms)
f2 <- psm(Surv(OS,event) ~ riskScore+PRS_type+Histology+Grade+Chemo_status+X1p19q_codeletion_status, data =  tb1, dist='lognormal')
med <- Quantile(f2)

surv <- Survival(f2)
nom <- nomogram(f2, fun=list(function(x) surv(1, x),
                             function(x) surv(3, x),
                             function(x) surv(5, x)),
                funlabel=c("1-year Survival Probability",
                           "3-year Survival Probability",
                           "5-year Survival Probability"))
plot(nom, xfrac=.6)

f2 <- psm(Surv(OS,event) ~ riskScore+PRS_type+Histology+Grade+Chemo_status+X1p19q_codeletion_status, data =  tb1,x = T,y = T,dist='lognormal')

cal1 <- calibrate(f2, 
                  cmethod='KM', 
                  method="boot", 
                  u=1, 
                  m=76, 
                  B=1000) 

cal1 <- calibrate(f2, 
                  cmethod = 'KM', 
                  method = "boot", 
                  u = 1, 
                  m = 200, 
                  B = 1000)

cal3 <- calibrate(f2, 
                  cmethod = 'KM', 
                  method = "boot", 
                  u = 3, 
                  m = 200, 
                  B = 1000)

cal5 <- calibrate(f2, 
                  cmethod = 'KM', 
                  method = "boot", 
                  u = 5, 
                  m = 200, 
                  B = 1000)

legend("bottomright", legend = c("1-Year", "3-Year", "5-Year"), col = c("red", "darkgreen", "darkblue"), lty = 1:3, lwd = 2)
plot(cal3,lwd=2,lty=1,
     conf.int=T,
     errbar.col="blue",
     col="blue",
     xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability of 1-Year DFS",
     ylab="Actual 1-Year DFS (proportion)",
     subtitles = T)
plot(cal1, lwd = 2, lty = 1, errbar.col = "red", col = "red", 
     xlim = c(0, 1), ylim = c(0, 1), xlab = "Nomogram-Predicted OS(%)", 
     ylab = "Actual OS(%)", subtitles = FALSE)

plot(cal3, lwd = 2, lty = 2, errbar.col = "green", col = "darkgreen", add = TRUE)
plot(cal5, lwd = 2, lty = 3, errbar.col = "blue", col = "darkblue", add = TRUE)

legend("bottomright", legend = c("1-Year", "3-Year", "5-Year"), col = c("red", "darkgreen", "darkblue"), lty = 1:3, lwd = 2)