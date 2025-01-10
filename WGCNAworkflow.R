rm(list = ls())
library(tibble)
library(dplyr)
library(tidyverse)
library(WGCNA)

enableWGCNAThreads()

options(stringsAsFactors = FALSE)

data<-read.csv("combined_4290_15824_108474.csv",header = T,row.names = 1)
pdata<-read.csv("20240522pdata.csv",header = T)
data<-as.data.frame(t(data))
data<-rownames_to_column(data,var = "sample")
data1<-inner_join(pdata,data,by="sample")
data1<-data1[,-c(2:5)]
rownames(data1)<-data1[,1]
data1<-data1[,-1]
data1<-as.data.frame(t(data1))
rownames(pdata)<-pdata[,1]
pdata<-pdata[,-1]

data1 <- as.data.frame(t(data1))

gsg = goodSamplesGenes(data1, verbose = 3)
if (!gsg$allOK) {
  
  if (sum(!gsg$goodGenes) > 0) 
    cat("Removing genes:", paste(colnames(data1)[!gsg$goodGenes], collapse = " "), "\n")
  if (sum(!gsg$goodSamples) > 0) 
    cat("Removing samples:", paste(rownames(data1)[!gsg$goodSamples], collapse = " "), "\n")
  
  data1 = data1[gsg$goodSamples, gsg$goodGenes]
}

threshold <- 1
data1 <- data1[, colMeans(data1) > threshold]

variance = apply(data1, 2, var)
cutoff = quantile(variance, 0.8) 
data1 = data1[, variance > cutoff]

cat("Filtered data contains", nrow(data1), "samples and", ncol(data1), "genes.\n")

write.table(data1, "filtered_data1.txt", sep = "\t", quote = FALSE, col.names = NA)

sampleTree = hclust(dist(data1), method = "average")
pdf(file = "1.sampleClustering.pdf", width = 15, height = 8)
par(cex = 0.6)
par(mar = c(0, 6, 0, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 2, cex.axis = 1.5, cex.main = 2)

dev.off()

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(data1, powerVector = powers, verbose = 5)

par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.90, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")

softPower = 9

adjacency = adjacency(data1, power = softPower)

TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

geneTree = hclust(as.dist(dissTOM), method = "average")

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
dynamicColors = labels2colors(dynamicMods)

sizeGrWindow(8, 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

MEList = moduleEigengenes(data1, colors = dynamicColors)
MEs = MEList$eigengenes

MEDiss = 1 - cor(MEs)

METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

mergeCutHeight = 0.25

merge = mergeCloseModules(data1, dynamicColors, cutHeight = mergeCutHeight, verbose = 3)

mergedColors = merge$colors
mergedMEs = merge$newMEs

moduleColors = mergedColors

sizeGrWindow(12, 9)
par(mar = c(5, 4, 4, 2) + 0.1)
plotDendroAndColors(geneTree, cbind(dynamicColors, moduleColors),
                    c("Dynamic Tree Cut", "Merged Modules"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs

samples <- rownames(data1)
traitData <- pdata[samples, ]

stopifnot(all(rownames(traitData) == rownames(data1)))

sampleDissimilarity = dist(data1)

sampleTree = hclust(sampleDissimilarity, method = "average")

traitColors = numbers2colors(traitData, signed = FALSE)

sizeGrWindow(12, 9)
par(mfrow = c(1, 1))
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(traitData),
                    main = "Sample dendrogram and trait heatmap",
                    dendroLabels = NULL, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

nGenes = ncol(data1)
nSamples = nrow(data1)
moduleTraitCor = cor(MEs, traitData, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

pdf(file="8_Module-trait relationships.pdf", width=10, height=10)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", 
                   signif(moduleTraitPvalue, 1), ")", sep = "") 
dim(textMatrix) = dim(moduleTraitCor) 
aaaa <- as.data.frame(traitData)
par(mar = c(6, 8.5, 3, 3)) 

labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels = names(aaaa), 
               yLabels = names(MEs), 
               ySymbols = names(MEs), 
               colorLabels = FALSE, 
               colors = greenWhiteRed(50), 
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.5, 
               zlim = c(-1, 1), 
               main = paste("Module-trait relationships"))
dev.off()

nSelect = 400
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]

selectTree = hclust(as.dist(selectTOM), method = "average")

selectColors = moduleColors[select]

plotDiss = selectTOM^7
diag(plotDiss) = NA

install.packages("gplots")
library(gplots)
pdf(file="Network_heatmap_plot_selected_genes.pdf", width = 9, height = 9)
mycol = colorpanel(250, "red", "orange", "lemonchiffon")
TOMplot(plotDiss, selectTree, selectColors, col = mycol, main = "Network heatmap plot selected genes")
dev.off()

pdf(file="Eigengene_dendrogram_2.pdf", width = 6, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0, 4, 2, 0), plotHeatmaps = FALSE)
dev.off()

pdf(file="Eigengene_adjacency_heatmap_2.pdf", width = 6, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3, 4, 2, 2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

modNames <- substring(names(MEs), 3)
modNames
geneModuleMembership <- as.data.frame(cor(data1, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")
geneModuleMembership[1:5, 1:5]

geneTraitSignificance <- as.data.frame(cor(data1, pdata, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(pdata), sep = "")
names(GSPvalue) <- paste("p.GS.", names(pdata), sep = "")
head(geneTraitSignificance)

module = "purple"
pheno = "GBM"
modNames = substring(names(MEs), 3)

module_column = match(module, modNames)
pheno_column = match(pheno, colnames(pdata))

moduleGenes <- moduleColors == module

par(mar = c(1, 3, 3, 3))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for LRG"),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(data1, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(data1)))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(data1, pdata, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(data1)))
names(geneTraitSignificance) <- paste("GS.", names(pdata), sep = "")
names(GSPvalue) <- paste("p.GS.", names(pdata), sep = "")

allModulesPValues <- list()

allModules <- unique(moduleColors)

for (module in allModules) {
  moduleGenes <- moduleColors == module
  moduleGeneNames <- colnames(data1)[moduleGenes]
  
  MM <- geneModuleMembership[moduleGenes, paste0("MM", module)]
  MMP <- MMPvalue[moduleGenes, paste0("p.MM", module)]
  
  GS <- geneTraitSignificance[moduleGenes, ]
  GSP <- GSPvalue[moduleGenes, ]
  
  moduleGeneInfo <- data.frame(Gene = moduleGeneNames, MM = MM, p.MM = MMP, GS = GS, p.GS = GSP)
  
  allModulesPValues[[module]] <- moduleGeneInfo
}

for (module in names(allModulesPValues)) {
  cat("Module:", module, "\n")
  print(head(allModulesPValues[[module]]))
  cat("\n")
}

for (module in names(allModulesPValues)) {
  write.table(allModulesPValues[[module]], file = paste0(module, "_module_gene_info.txt"), sep = "\t", row.names = FALSE)
}
TOM = TOMsimilarityFromExpr(data1, power = 9)
module = "brown"
probes = colnames(data1) 
inModule = (moduleColors==module)
modProbes = probes[inModule]
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.1,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
)

library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(tibble)
green<-read.csv("hubgene_MMGS_green.csv")

green1 <- bitr(green$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene<-green1$ENTREZID
ego_ALL <- enrichGO(gene = gene,
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    minGSSize = 1,
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05,
                    readable = TRUE)
ego_ALL<-as.data.frame(ego_ALL@result)
