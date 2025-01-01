rm(list = ls())
## Sample reduction
data<-read.csv("combined_4290_15824_108474.csv",header = T,row.names = 1)
pdata<-read.csv("20240522pdata.csv",header = T)
data<-as.data.frame(t(data))
library(tibble)
data<-rownames_to_column(data,var = "sample")
library(dplyr)
library(tidyverse)
data1<-inner_join(pdata,data,by="sample")
data1<-data1[,-c(2:5)]
rownames(data1)<-data1[,1]
data1<-data1[,-1]
data1<-as.data.frame(t(data1))
rownames(pdata)<-pdata[,1]
pdata<-pdata[,-1]
enableWGCNAThreads()

#load the packages
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  install.packages("WGCNA")
}

library(WGCNA)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("GO.db")
BiocManager::install("impute")

# Set environment parameters
options(stringsAsFactors = FALSE)

# Transpose the expression matrix so that rows are samples and columns are genes
data1 <- as.data.frame(t(data1))


# Check data quality
gsg = goodSamplesGenes(data1, verbose = 3)
if (!gsg$allOK) {
  # If there are genes or samples that do not meet the standards, remove them
  if (sum(!gsg$goodGenes) > 0) 
    cat("Removing genes:", paste(colnames(data1)[!gsg$goodGenes], collapse = " "), "\n")
  if (sum(!gsg$goodSamples) > 0) 
    cat("Removing samples:", paste(rownames(data1)[!gsg$goodSamples], collapse = " "), "\n")
  
  data1 = data1[gsg$goodSamples, gsg$goodGenes]
}

# Filter low-expression genes
threshold <- 1
data1 <- data1[, colMeans(data1) > threshold]

# Filter low-variance genes
variance = apply(data1, 2, var)
cutoff = quantile(variance, 0.8) # Adjust as needed
data1 = data1[, variance > cutoff]

# Print the number of genes and samples remaining after filtering
cat("Filtered data contains", nrow(data1), "samples and", ncol(data1), "genes.\n")

# Save the filtered data
write.table(data1, "filtered_data1.txt", sep = "\t", quote = FALSE, col.names = NA)

# Sample selection, remove outlier samples
sampleTree = hclust(dist(data1), method = "average")
pdf(file = "1.sampleClustering.pdf", width = 15, height = 8)
par(cex = 0.6)
par(mar = c(0, 6, 0, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", 
     xlab = "", cex.lab = 2, cex.axis = 1.5, cex.main = 2)
# abline(h = 180, col = "red") # Adjust the cutting height as needed
dev.off()

# Select the appropriate soft threshold (power)
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(data1, powerVector = powers, verbose = 5)

# Set the graphics window size
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# Plot the relationship between the R² of the scale-free topology fitting index and different power values
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.90, col = "red") # If R² reaches 0.90, it is considered a scale-free network

# Plot the relationship between average connectivity and different power values
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     labels = powers, cex = 0.9, col = "red")

# Select the optimal soft threshold (e.g., power = 6)
softPower = 9

# Calculate the adjacency matrix
adjacency = adjacency(data1, power = softPower)

# Calculate the topological overlap matrix (TOM)
TOM = TOMsimilarity(adjacency)
dissTOM = 1 - TOM

# Hierarchical clustering
geneTree = hclust(as.dist(dissTOM), method = "average")

# Dynamic tree cutting
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, 
                            deepSplit = 2, pamRespectsDendro = FALSE, 
                            minClusterSize = 30)
dynamicColors = labels2colors(dynamicMods)

# Plot the clustering tree and module colors
sizeGrWindow(8, 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Calculate the module eigengenes
MEList = moduleEigengenes(data1, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate the correlation between module eigengenes
MEDiss = 1 - cor(MEs)

# Plot the dendrogram of module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")

# Set the merge threshold (e.g., 0.25)
mergeCutHeight = 0.25

# Merge similar modules
merge = mergeCloseModules(data1, dynamicColors, cutHeight = mergeCutHeight, verbose = 3)

# Get the module eigengenes after merging
mergedColors = merge$colors
mergedMEs = merge$newMEs

# Update module colors
moduleColors = mergedColors

# Plot the merged clustering tree and module colors
sizeGrWindow(12, 9)
par(mar = c(5, 4, 4, 2) + 0.1)
plotDendroAndColors(geneTree, cbind(dynamicColors, moduleColors),
                    c("Dynamic Tree Cut", "Merged Modules"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Set module colors and module eigengenes
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs

# Read trait data pdata
# Ensure the sample order is consistent
samples <- rownames(data1)
traitData <- pdata[samples, ]

# Check if the trait data and expression data match
stopifnot(all(rownames(traitData) == rownames(data1)))

# Calculate the sample distance matrix
sampleDissimilarity = dist(data1)

# Perform hierarchical clustering
sampleTree = hclust(sampleDissimilarity, method = "average")

# Convert trait data into a color matrix
traitColors = numbers2colors(traitData, signed = FALSE)

# Plot the sample dendrogram and trait heatmap
sizeGrWindow(12, 9)
par(mfrow = c(1, 1))
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(traitData),
                    main = "Sample dendrogram and trait heatmap",
                    dendroLabels = NULL, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Calculate the correlation between modules and traits
nGenes = ncol(data1)
nSamples = nrow(data1)
moduleTraitCor = cor(MEs, traitData, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Plot the heatmap of module-trait correlations
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

# Generate TOM plot
nSelect = 400
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]

# Gene clustering
selectTree = hclust(as.dist(selectTOM), method = "average")

# Get the module colors for the selected genes
selectColors = moduleColors[select]

# Set TOM plot parameters
plotDiss = selectTOM^7
diag(plotDiss) = NA

# Plot the TOM plot
install.packages("gplots")
library(gplots)
pdf(file="Network_heatmap_plot_selected_genes.pdf", width = 9, height = 9)
mycol = colorpanel(250, "red", "orange", "lemonchiffon")
TOMplot(plotDiss, selectTree, selectColors, col = mycol, 
        main = "Network heatmap plot selected genes")
dev.off()

# Plot the dendrogram and heatmap of module eigengenes
pdf(file="Eigengene_dendrogram_2.pdf", width = 6, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", 
                      marDendro = c(0, 4, 2, 0), plotHeatmaps = FALSE)
dev.off()

pdf(file="Eigengene_adjacency_heatmap_2.pdf", width = 6, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", 
                      marHeatmap = c(3, 4, 2, 2), plotDendrograms = FALSE, 
                      xLabelsAngle = 90)
dev.off()

# Module gene extraction
modNames <- substring(names(MEs), 3)
modNames
geneModuleMembership <- as.data.frame(cor(data1, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                           nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")
geneModuleMembership[1:5, 1:5]

# Calculate the gene-trait correlation matrix 
# Only continuous traits can be calculated; if discrete variables are present, they should be converted to 0-1 matrices when constructing the sample table.
geneTraitSignificance <- as.data.frame(cor(data1, pdata, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                           nSamples))
names(geneTraitSignificance) <- paste("GS.", names(pdata), sep = "")
names(GSPvalue) <- paste("p.GS.", names(pdata), sep = "")
head(geneTraitSignificance)

# Finally, combine the two correlation matrices, specify the module of interest for analysis
module = "purple"
pheno = "GBM"
modNames = substring(names(MEs), 3)

# Get the columns of interest
module_column = match(module, modNames)
pheno_column = match(pheno, colnames(pdata))

# Get the genes within the module
moduleGenes <- moduleColors == module

# Visualize the correlation between genes and modules, phenotypes, and plot scatter plots
par(mar = c(1, 3, 3, 3))

verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for LRG"),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# All module gene extraction
modNames <- substring(names(MEs), 3)

# Calculate gene module membership (MM)
geneModuleMembership <- as.data.frame(cor(data1, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                           nrow(data1)))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

# Calculate gene trait correlation (GS)
geneTraitSignificance <- as.data.frame(cor(data1, pdata, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                           nrow(data1)))
names(geneTraitSignificance) <- paste("GS.", names(pdata), sep = "")
names(GSPvalue) <- paste("p.GS.", names(pdata), sep = "")

# Initialize a list to store the results
allModulesPValues <- list()

# Get all module colors
allModules <- unique(moduleColors)

# Iterate through all modules and extract genes and associated p-values
for (module in allModules) {
  # Get the genes within the module
  moduleGenes <- moduleColors == module
  moduleGeneNames <- colnames(data1)[moduleGenes]
  # Get module membership and p-values
  MM <- geneModuleMembership[moduleGenes, paste0("MM", module)]
  MMP <- MMPvalue[moduleGenes, paste0("p.MM", module)]
  # Get gene trait correlation and p-values
  GS <- geneTraitSignificance[moduleGenes, ]
  GSP <- GSPvalue[moduleGenes, ]
  # Create a data frame to store gene names, MM, MM p-values, GS, and GS p-values
  moduleGeneInfo <- data.frame(Gene = moduleGeneNames, MM = MM, 
                               p.MM = MMP, GS = GS, p.GS = GSP)
  # Store the data frame in the list
  allModulesPValues[[module]] <- moduleGeneInfo
}

# Print the genes and associated p-values for all modules
for (module in names(allModulesPValues)) {
  cat("Module:", module, "\n")
  print(head(allModulesPValues[[module]]))
  cat("\n")
}

# If you need to save the results to a file
for (module in names(allModulesPValues)) {
  write.table(allModulesPValues[[module]], 
              file = paste0(module, "_module_gene_info.txt"), 
              sep = "\t", row.names = FALSE)
}
#####Export genes to cytoscape
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(data1, power = 9); 
# Select module
module = "brown"
# Select module probes
probes = colnames(data1) ## In our example, the probe is the gene
inModule = (moduleColors==module)
modProbes = probes[inModule]
## Also extract the names of the genes in the specified module
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
## Gene relationship matrix corresponding to the module
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), 
                   ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), 
                   ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.1,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
)

####Hub gene GO analysis
library(AnnotationDbi)
library(org.Hs.eg.db)#Gene annotation package
library(clusterProfiler)#Enrichment package
library(dplyr)
library(ggplot2)#Plotting package
library(tibble)
green<-read.csv("hubgene_MMGS_green.csv")

# Gene ID conversion, the ID type used for KEGG and GO enrichment is ENTREZID）
green1 <- bitr(green$X, fromType = "SYMBOL", toType = "ENTREZID", 
               OrgDb = org.Hs.eg.db)
gene<-green1$ENTREZID
ego_ALL <- enrichGO(gene = gene,#We defined it above
                    OrgDb=org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",#Types of GO enrichment
                    pAdjustMethod = "BH",#This doesn't need to be managed, BH is generally used
                    minGSSize = 1,
                    pvalueCutoff = 0.01,#P value can be 0.05
                    qvalueCutoff = 0.05,
                    readable = TRUE)
ego_ALL<-as.data.frame(ego_ALL@result)

