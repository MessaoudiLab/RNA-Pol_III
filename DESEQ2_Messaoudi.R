#!/bin/R

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)

## ----results='hide', message=FALSE-----------------------------------------------------------
# load the libraries
print(.libPaths())
library(DESeq2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(BiocParallel)
#library(org.Hs.eg.db)


## --------------------------------------------------------------------------------------------
# Read in the raw read counts
rawCounts <- read.delim(args[1])

  
# Read in the sample mappings
sampleData <- read.csv(args[2])

# Retrieve GeneIDs
# Retrieved April 2022
#gns <- select(org.Hs.eg.db, rawCounts$Geneid, "SYMBOL", "ENSEMBL")
gns <- read.table(args[3])
names(gns) <- c("Geneid", "symbol")



## --------------------------------------------------------------------------------------------
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
  
}

justCounts <- rawCounts[7:ncol(rawCounts)]
justCounts[is.na(justCounts)] <- 0
genes <- rawCounts %>%
  dplyr::select(Geneid, Length)
genes[is.na(genes)] <- 1
tpms <- apply(justCounts, 2, function(x) tpm(x, genes$Length))
tpms <- as.data.frame(tpms)
tpms$Geneid <- genes$Geneid

tpms.annotated <- merge(gns, tpms)



## --------------------------------------------------------------------------------------------
#Filter for only untransfected cell lines
sampleData.untransfected <- sampleData %>%
  filter(Transfection == 'Untransfected')


## --------------------------------------------------------------------------------------------
# Filter for A549 (WT) cell line
# Remove nR115 samples to avoid batch effect
sampleData.untransfected.wt.nr115Removed <- sampleData.untransfected %>%
  filter(Cell_line == 'A549') %>%
  filter(Batch=='nR052')

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.untransfected.wt <- sampleData.untransfected.wt.nr115Removed$Filename
rawCounts.matrix.untransfected.wt <- as.matrix(rawCounts[,sampleIndex.untransfected.wt])
rownames(rawCounts.matrix.untransfected.wt) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.untransfected.wt.nr115Removed) <- sampleData.untransfected.wt.nr115Removed$Filename
keep <- c("Condition", "SampleName")
sampleData.untransfected.wt.nr115Removed <- sampleData.untransfected.wt.nr115Removed[,keep]
sampleData.untransfected.wt.nr115Removed$Condition <- factor(sampleData.untransfected.wt.nr115Removed$Condition)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.untransfected.wt.ordered <- rawCounts.matrix.untransfected.wt[,unique(rownames(sampleData.untransfected.wt.nr115Removed))]
all(colnames(rawCounts.matrix.untransfected.wt.ordered) == rownames(sampleData.untransfected.wt.nr115Removed))

# Set the reference level for each factor
sampleData.untransfected.wt.nr115Removed$Condition <- relevel(sampleData.untransfected.wt.nr115Removed$Condition, ref='NonInfected')

#Set NA values to 0
rawCounts.matrix.untransfected.wt.ordered[is.na(rawCounts.matrix.untransfected.wt.ordered)] <- 0


## --------------------------------------------------------------------------------------------
# Create the DEseq2DataSet object for DE analysis
deseq2Data.untransfected.wt <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.untransfected.wt.ordered, colData=sampleData.untransfected.wt.nr115Removed, design= ~Condition)


## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.untransfected.wt)
plotPCA(rld, intgroup=c('Condition'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.untransfected.wt)
dim(deseq2Data.untransfected.wt[rowSums(counts(deseq2Data.untransfected.wt)) > 5, ])
deseq2Data.untransfected.wt <- deseq2Data.untransfected.wt[rowSums(counts(deseq2Data.untransfected.wt)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.untransfected.wt <- DESeq(deseq2Data.untransfected.wt, parallel = TRUE)

# Extract differential expression results
deseq2Results.untransfected.wt <- results(deseq2Data.untransfected.wt, alpha=0.05, contrast=c("Condition", "Infected", "NonInfected"))
deseq2Results.untransfected.wt.tidy <- results(deseq2Data.untransfected.wt, alpha=0.05, contrast=c("Condition", "Infected", "NonInfected"), tidy=TRUE)

# View summary of results
summary(deseq2Results.untransfected.wt)

gns.deseq2results.untransfected.wt <- gns
names(gns.deseq2results.untransfected.wt) <- c("row", "symbol")
df_results.untransfected.wt <- merge(deseq2Results.untransfected.wt.tidy, gns.deseq2results.untransfected.wt, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.untransfected.wt <- which(deseq2Results.untransfected.wt$padj < 0.05 & abs(deseq2Results.untransfected.wt$log2FoldChange) > 1.1) 

# Get the annotation column labels
df.untransfected.wt <- as.data.frame(colData(deseq2Data.untransfected.wt)[,c("Condition")])

# Perform the VST transformation of the counts
vst.untransfected.wt <- vst(deseq2Data.untransfected.wt)

# Match the rownames and column names
rownames(df.untransfected.wt) <- colnames(vst.untransfected.wt)
colnames(df.untransfected.wt) <- c("Condition")

gns.untransfected.wt <- gns %>%
  filter(Geneid %in% row.names(vst.untransfected.wt))
row.names(vst.untransfected.wt)[match(gns.untransfected.wt[,1], row.names(vst.untransfected.wt))] <- gns.untransfected.wt[,2]

# Create Heat Map
pheatmap(assay(vst.untransfected.wt)[select.untransfected.wt,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df.untransfected.wt )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method, plot MA 
plotMA(deseq2Results.untransfected.wt)


## --------------------------------------------------------------------------------------------
# Filter for MAVS-KO cell line
sampleData.untransfected.mavs <- sampleData.untransfected %>%
  filter(Cell_line == 'MAVS-KO')

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.untransfected.mavs <- sampleData.untransfected.mavs$Filename
rawCounts.matrix.untransfected.mavs <- as.matrix(rawCounts[,sampleIndex.untransfected.mavs])
rownames(rawCounts.matrix.untransfected.mavs) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.untransfected.mavs) <- sampleData.untransfected.mavs$Filename
keep <- c("Condition", "SampleName")
sampleData.untransfected.mavs <- sampleData.untransfected.mavs[,keep]
sampleData.untransfected.mavs$Condition <- factor(sampleData.untransfected.mavs$Condition)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.untransfected.mavs.ordered <- rawCounts.matrix.untransfected.mavs[,unique(rownames(sampleData.untransfected.mavs))]
all(colnames(rawCounts.matrix.untransfected.mavs.ordered) == rownames(sampleData.untransfected.mavs))

# Set the reference level for each factor
sampleData.untransfected.mavs$Condition <- relevel(sampleData.untransfected.mavs$Condition, ref='NonInfected')

#Set NA values to 0
rawCounts.matrix.untransfected.mavs.ordered[is.na(rawCounts.matrix.untransfected.mavs.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.untransfected.mavs <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.untransfected.mavs.ordered, colData=sampleData.untransfected.mavs, design= ~ Condition)


## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.untransfected.mavs)
plotPCA(rld, intgroup=c('Condition'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.untransfected.mavs)
dim(deseq2Data.untransfected.mavs[rowSums(counts(deseq2Data.untransfected.mavs)) > 5, ])
deseq2Data.untransfected.mavs <- deseq2Data.untransfected.mavs[rowSums(counts(deseq2Data.untransfected.mavs)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.untransfected.mavs <- DESeq(deseq2Data.untransfected.mavs, parallel = TRUE)

# Extract differential expression results
deseq2Results.untransfected.mavs <- results(deseq2Data.untransfected.mavs, alpha=0.05, contrast=c("Condition", "Infected", "NonInfected"))
deseq2Results.untransfected.mavs.tidy <- results(deseq2Data.untransfected.mavs, alpha=0.05, contrast=c("Condition", "Infected", "NonInfected"), tidy=TRUE)

# View summary of results
summary(deseq2Results.untransfected.mavs)

# Set GeneIDs
gns.deseq2results.untransfected.mavs <- gns 
names(gns.deseq2results.untransfected.mavs) <- c("row", "symbol")
df_results.untransfected.mavs <- merge(deseq2Results.untransfected.mavs.tidy, gns.deseq2results.untransfected.mavs, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.untransfected.mavs <- which(deseq2Results.untransfected.mavs$padj < 0.001 & abs(deseq2Results.untransfected.mavs$log2FoldChange) > 3) 

# Get the annotation column labels
df.untransfected.mavs <- as.data.frame(colData(deseq2Data.untransfected.mavs)[,c("Condition")])

# Perform the VST transformation of the counts
vst.untransfected.mavs <- vst(deseq2Data.untransfected.mavs)

# Match the rownames and column names
rownames(df.untransfected.mavs) <- colnames(vst.untransfected.mavs)
colnames(df.untransfected.mavs) <- c("Condition")

gns.untransfected.mavs <- gns %>%
  filter(Geneid %in% row.names(vst.untransfected.mavs))
row.names(vst.untransfected.mavs)[match(gns.untransfected.mavs[,1], row.names(vst.untransfected.mavs))] <- gns.untransfected.mavs[,2]

# Create Heat Map
pheatmap(assay(vst.untransfected.mavs)[select.untransfected.mavs,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df.untransfected.mavs )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.untransfected.mavs)


## --------------------------------------------------------------------------------------------
# Filter for PKR-KO cell line
sampleData.untransfected.pkr <- sampleData.untransfected %>%
  filter(Cell_line == 'PKR-KO')

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.untransfected.pkr <- sampleData.untransfected.pkr$Filename
rawCounts.matrix.untransfected.pkr <- as.matrix(rawCounts[,sampleIndex.untransfected.pkr])
rownames(rawCounts.matrix.untransfected.pkr) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.untransfected.pkr) <- sampleData.untransfected.pkr$Filename
keep <- c("Condition", "SampleName")
sampleData.untransfected.pkr <- sampleData.untransfected.pkr[,keep]
sampleData.untransfected.pkr$Condition <- factor(sampleData.untransfected.pkr$Condition)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.untransfected.pkr.ordered <- rawCounts.matrix.untransfected.pkr[,unique(rownames(sampleData.untransfected.pkr))]
all(colnames(rawCounts.matrix.untransfected.pkr.ordered) == rownames(sampleData.untransfected.pkr))

# Set the reference level for each factor
sampleData.untransfected.pkr$Condition <- relevel(sampleData.untransfected.pkr$Condition, ref='NonInfected')

#Set NA values to 0
rawCounts.matrix.untransfected.pkr.ordered[is.na(rawCounts.matrix.untransfected.pkr.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.untransfected.pkr <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.untransfected.pkr.ordered, colData=sampleData.untransfected.pkr, design= ~ Condition)


## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.untransfected.pkr)
plotPCA(rld, intgroup=c('Condition'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.untransfected.pkr)
dim(deseq2Data.untransfected.pkr[rowSums(counts(deseq2Data.untransfected.pkr)) > 5, ])
deseq2Data.untransfected.pkr <- deseq2Data.untransfected.pkr[rowSums(counts(deseq2Data.untransfected.pkr)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.untransfected.pkr <- DESeq(deseq2Data.untransfected.pkr, parallel = TRUE)

# Extract differential expression results
deseq2Results.untransfected.pkr <- results(deseq2Data.untransfected.pkr, alpha=0.05, contrast=c("Condition", "Infected", "NonInfected"))
deseq2Results.untransfected.pkr.tidy <- results(deseq2Data.untransfected.pkr, alpha=0.05, contrast=c("Condition", "Infected", "NonInfected"), tidy=TRUE)

# View summary of results
summary(deseq2Results.untransfected.pkr)


# Retrieve GeneIDs
gns.deseq2results.untransfected.pkr <- gns
names(gns.deseq2results.untransfected.pkr) <- c("row", "symbol")
df_results.untransfected.pkr <- merge(deseq2Results.untransfected.pkr.tidy, gns.deseq2results.untransfected.pkr, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.untransfected.pkr <- which(deseq2Results.untransfected.pkr$padj < 0.05 & abs(deseq2Results.untransfected.pkr$log2FoldChange)>4) 
# Get the annotation column labels
df.untransfected.pkr <- as.data.frame(colData(deseq2Data.untransfected.pkr)[,c("Condition")])
# Perform the VST transformation of the counts
vst.untransfected.pkr <- vst(deseq2Data.untransfected.pkr)
# Match the row names and column names
rownames(df.untransfected.pkr) <- colnames(vst.untransfected.pkr)
colnames(df.untransfected.pkr) <- c("Condition")

# get Geneids
gns.untransfected.pkr <- gns %>%
  filter(Geneid %in% row.names(vst.untransfected.pkr))
row.names(vst.untransfected.pkr)[match(gns.untransfected.pkr[,1], row.names(vst.untransfected.pkr))] <- gns.untransfected.pkr[,2]

# Create Heat Map
pheatmap(assay(vst.untransfected.pkr)[select.untransfected.pkr,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.untransfected.pkr )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.untransfected.pkr)


## --------------------------------------------------------------------------------------------
# Filter for MAVS-KO cell line
sampleData.untransfected.rl <- sampleData.untransfected %>%
  filter(Cell_line == 'RL-KO')

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.untransfected.rl <- sampleData.untransfected.rl$Filename
rawCounts.matrix.untransfected.rl <- as.matrix(rawCounts[,sampleIndex.untransfected.rl])
rownames(rawCounts.matrix.untransfected.rl) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.untransfected.rl) <- sampleData.untransfected.rl$Filename
keep <- c("Condition", "SampleName")
sampleData.untransfected.rl <- sampleData.untransfected.rl[,keep]
sampleData.untransfected.rl$Condition <- factor(sampleData.untransfected.rl$Condition)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.untransfected.rl.ordered <- rawCounts.matrix.untransfected.rl[,unique(rownames(sampleData.untransfected.rl))]
all(colnames(rawCounts.matrix.untransfected.rl.ordered) == rownames(sampleData.untransfected.rl))

# Set the reference level for each factor
sampleData.untransfected.rl$Condition <- relevel(sampleData.untransfected.rl$Condition, ref='NonInfected')

#Set NA values to 0
rawCounts.matrix.untransfected.rl.ordered[is.na(rawCounts.matrix.untransfected.rl.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.untransfected.rl <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.untransfected.rl.ordered, colData=sampleData.untransfected.rl, design= ~ Condition)


## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.untransfected.rl)
plotPCA(rld, intgroup=c('Condition'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.untransfected.rl)
dim(deseq2Data.untransfected.rl[rowSums(counts(deseq2Data.untransfected.rl)) > 5, ])
deseq2Data.untransfected.rl <- deseq2Data.untransfected.rl[rowSums(counts(deseq2Data.untransfected.rl)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.untransfected.rl <- DESeq(deseq2Data.untransfected.rl, parallel = TRUE)

# Extract differential expression results
deseq2Results.untransfected.rl <- results(deseq2Data.untransfected.rl, alpha=0.05, contrast=c("Condition", "Infected", "NonInfected"))
deseq2Results.untransfected.rl.tidy <- results(deseq2Data.untransfected.rl, alpha=0.05, contrast=c("Condition", "Infected", "NonInfected"), tidy=TRUE)

# View summary of results
summary(deseq2Results.untransfected.rl)

# Get geneIDs
gns.deseq2results.untransfected.rl <- gns
names(gns.deseq2results.untransfected.rl) <- c("row", "symbol")
df_results.untransfected.rl <- merge(deseq2Results.untransfected.rl.tidy, gns.deseq2results.untransfected.rl, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.untransfected.rl <- which(deseq2Results.untransfected.rl$padj < 0.05 & abs(deseq2Results.untransfected.rl$log2FoldChange) > 4)

# Get the annotation column labels
df.untransfected.rl <- as.data.frame(colData(deseq2Data.untransfected.rl)[,c("Condition")])

# Perform the VST transformation of the counts
vst.untransfected.rl <- vst(deseq2Data.untransfected.rl)

# Match the rownames and column names
rownames(df.untransfected.rl) <- colnames(vst.untransfected.rl)
colnames(df.untransfected.rl) <- c("Condition")

# get geneids
gns.untransfected.rl <- gns %>%
  filter(Geneid %in% row.names(vst.untransfected.rl))
row.names(vst.untransfected.rl)[match(gns.untransfected.rl[,1], row.names(vst.untransfected.rl))] <- gns.untransfected.rl[,2]

# Create Heat Map
pheatmap(assay(vst.untransfected.rl)[select.untransfected.rl,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.untransfected.rl )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.untransfected.rl)


## --------------------------------------------------------------------------------------------
# Filter for TLR3-KO cell line
sampleData.untransfected.tlr3 <- sampleData.untransfected %>%
  filter(Cell_line == 'TLR3-KO')

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex <- sampleData.untransfected.tlr3$Filename
rawCounts.matrix.untransfected.tlr3 <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts.matrix.untransfected.tlr3) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.untransfected.tlr3) <- sampleData.untransfected.tlr3$Filename
keep <- c("Condition", "SampleName")
sampleData.untransfected.tlr3 <- sampleData.untransfected.tlr3[,keep]
sampleData.untransfected.tlr3$Condition <- factor(sampleData.untransfected.tlr3$Condition)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.untransfected.tlr3.ordered <- rawCounts.matrix.untransfected.tlr3[,unique(rownames(sampleData.untransfected.tlr3))]
all(colnames(rawCounts.matrix.untransfected.tlr3.ordered) == rownames(sampleData.untransfected.tlr3))

# Set the reference level for each factor
sampleData.untransfected.tlr3$Condition <- relevel(sampleData.untransfected.tlr3$Condition, ref='NonInfected')

#Set NA values to 0
rawCounts.matrix.untransfected.tlr3.ordered[is.na(rawCounts.matrix.untransfected.tlr3.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.untransfected.tlr3 <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.untransfected.tlr3.ordered, colData=sampleData.untransfected.tlr3, design= ~ Condition)


## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.untransfected.tlr3)
plotPCA(rld, intgroup=c('Condition'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.untransfected.tlr3)
dim(deseq2Data.untransfected.tlr3[rowSums(counts(deseq2Data.untransfected.tlr3)) > 5, ])
deseq2Data.untransfected.tlr3 <- deseq2Data.untransfected.tlr3[rowSums(counts(deseq2Data.untransfected.tlr3)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.untransfected.tlr3 <- DESeq(deseq2Data.untransfected.tlr3, parallel = TRUE)

# Extract differential expression results
deseq2Results.untransfected.tlr3 <- results(deseq2Data.untransfected.tlr3, alpha=0.05, contrast=c("Condition", "Infected", "NonInfected"))
deseq2Results.untransfected.tlr3.tidy <- results(deseq2Data.untransfected.tlr3, alpha=0.05, contrast=c("Condition", "Infected", "NonInfected"), tidy=TRUE)

# View summary of results
summary(deseq2Results.untransfected.tlr3)

# Get geneIDs 
gns.deseq2results.untransfected.tlr3 <- gns
names(gns.deseq2results.untransfected.tlr3) <- c("row", "symbol")
df_results.untransfected.tlr3 <- merge(deseq2Results.untransfected.tlr3.tidy, gns.deseq2results.untransfected.tlr3, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.untransfected.tlr3 <- which(deseq2Results.untransfected.tlr3$padj < 0.05 & abs(deseq2Results.untransfected.tlr3$log2FoldChange)>4)

# Get the annotation column labels
df.untransfected.tlr3 <- as.data.frame(colData(deseq2Data.untransfected.tlr3)[,c("Condition")])

# Perform the VST transformation of the counts
vst.untransfected.tlr3 <- vst(deseq2Data.untransfected.tlr3)

# Match the rownames and column names
rownames(df.untransfected.tlr3) <- colnames(vst.untransfected.tlr3)
colnames(df.untransfected.tlr3) <- c("Condition")

# get geneids
gns.untransfected.tlr3 <- gns %>%
  filter(Geneid %in% row.names(vst.untransfected.tlr3))
row.names(vst.untransfected.tlr3)[match(gns.untransfected.tlr3[,1], row.names(vst.untransfected.tlr3))] <- gns.untransfected.tlr3[,2]

# Create Heat Map
pheatmap(assay(vst.untransfected.tlr3)[select.untransfected.tlr3,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.untransfected.tlr3 )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.untransfected.tlr3)


## --------------------------------------------------------------------------------------------
#Filter for only transfected cells
sampleData.transfected <- sampleData %>%
  filter(Transfection != 'Untransfected')


## --------------------------------------------------------------------------------------------
# Filter for A549 (WT) cell line and well aligned samples
sampleData.transfected.wt <- sampleData.transfected %>%
  filter(Cell_line == 'A549') %>%
  filter(Batch == 'nR052')
  

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.transfected.wt <- sampleData.transfected.wt$Filename
rawCounts.matrix.transfected.wt <- as.matrix(rawCounts[,sampleIndex.transfected.wt])
rownames(rawCounts.matrix.transfected.wt) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.transfected.wt) <- sampleData.transfected.wt$Filename
keep <- c("Condition", "SampleName", "Batch", "Transfection")
sampleData.transfected.wt <- sampleData.transfected.wt[,keep]
sampleData.transfected.wt$Condition <- factor(sampleData.transfected.wt$Condition)
sampleData.transfected.wt$Transfection <- factor(sampleData.transfected.wt$Transfection)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.transfected.wt.ordered <- rawCounts.matrix.transfected.wt[,unique(rownames(sampleData.transfected.wt))]
all(colnames(rawCounts.matrix.transfected.wt.ordered) == rownames(sampleData.transfected.wt))

# Set the reference level for each factor
sampleData.transfected.wt$Condition <- relevel(sampleData.transfected.wt$Condition, ref='NonInfected')
sampleData.transfected.wt$Transfection <- relevel(sampleData.transfected.wt$Transfection, ref='siCtrl')

#Set NA values to 0
rawCounts.matrix.transfected.wt.ordered[is.na(rawCounts.matrix.transfected.wt.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.transfected.wt <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.transfected.wt.ordered, colData=sampleData.transfected.wt, design= ~Transfection*Condition)




## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.transfected.wt)
plotPCA(rld, intgroup=c('Condition', 'Transfection'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.transfected.wt)
dim(deseq2Data.transfected.wt[rowSums(counts(deseq2Data.transfected.wt)) > 5, ])
deseq2Data.transfected.wt <- deseq2Data.transfected.wt[rowSums(counts(deseq2Data.transfected.wt)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.transfected.wt <- DESeq(deseq2Data.transfected.wt, parallel = TRUE)

# Extract differential expression results
deseq2Results.transfected.wt <- results(deseq2Data.transfected.wt, alpha=0.05,list(c("Condition_Infected_vs_NonInfected", "TransfectionsiBrf1.ConditionInfected")))
deseq2Results.transfected.wt.tidy <- results(deseq2Data.transfected.wt, alpha=0.05,list(c("Condition_Infected_vs_NonInfected", "TransfectionsiBrf1.ConditionInfected")), tidy=TRUE)


# View summary of results
summary(deseq2Results.transfected.wt)

# Get geneids
gns.deseq2Results.transfected.wt <- gns
names(gns.deseq2Results.transfected.wt) <- c("row", "symbol")
df_results.transfected.wt <- merge(deseq2Results.transfected.wt.tidy, gns.deseq2Results.transfected.wt, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.transfected.wt <- which(deseq2Results.transfected.wt$padj < 0.00000000005) 
df.transfected.wt <- as.data.frame(colData(deseq2Data.transfected.wt)[,c("Condition", "Transfection")])
vst.transfected.wt <- vst(deseq2Data.transfected.wt)
rownames(df.transfected.wt) <- colnames(vst.transfected.wt)
colnames(df.transfected.wt) <- c("Condition", "Transfection")

# get geneids 
gns.transfected.wt <- gns %>%
  filter(Geneid %in% row.names(vst.transfected.wt))
row.names(vst.transfected.wt)[match(gns.transfected.wt[,1], row.names(vst.transfected.wt))] <- gns.transfected.wt[,2]

pheatmap(assay(vst.transfected.wt)[select.transfected.wt,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.transfected.wt )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.transfected.wt)


## --------------------------------------------------------------------------------------------
# Filter for MAVS-KO cell line
sampleData.transfected.mavs <- sampleData.transfected %>%
  filter(Cell_line == 'MAVS-KO') 


# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.transfected.mavs <- sampleData.transfected.mavs$Filename
rawCounts.matrix.transfected.mavs <- as.matrix(rawCounts[,sampleIndex.transfected.mavs])
rownames(rawCounts.matrix.transfected.mavs) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.transfected.mavs) <- sampleData.transfected.mavs$Filename
keep <- c("Condition", "SampleName","Transfection")
sampleData.transfected.mavs <- sampleData.transfected.mavs[,keep]
sampleData.transfected.mavs$Condition <- factor(sampleData.transfected.mavs$Condition)
sampleData.transfected.mavs$Transfection <- factor(sampleData.transfected.mavs$Transfection)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.transfected.mavs.ordered <- rawCounts.matrix.transfected.mavs[,unique(rownames(sampleData.transfected.mavs))]
all(colnames(rawCounts.matrix.transfected.mavs.ordered) == rownames(sampleData.transfected.mavs))

# Set the reference level for each factor
sampleData.transfected.mavs$Condition <- relevel(sampleData.transfected.mavs$Condition, ref='NonInfected')
sampleData.transfected.mavs$Transfection <- relevel(sampleData.transfected.mavs$Transfection, ref='siCtrl')

#Set NA values to 0
rawCounts.matrix.transfected.mavs.ordered[is.na(rawCounts.matrix.transfected.mavs.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.transfected.mavs <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.transfected.mavs.ordered, colData=sampleData.transfected.mavs, design= ~Transfection*Condition)




## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.transfected.mavs)
plotPCA(rld, intgroup=c('Condition', 'Transfection'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.transfected.mavs)
dim(deseq2Data.transfected.mavs[rowSums(counts(deseq2Data.transfected.mavs)) > 5, ])
deseq2Data.transfected.mavs <- deseq2Data.transfected.mavs[rowSums(counts(deseq2Data.transfected.mavs)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.transfected.mavs <- DESeq(deseq2Data.transfected.mavs, parallel = TRUE)

# Extract differential expression results
deseq2Results.transfected.mavs <- results(deseq2Data.transfected.mavs, alpha=0.05,list(c("Condition_Infected_vs_NonInfected", "TransfectionsiBrf1.ConditionInfected")))
deseq2Results.transfected.mavs.tidy  <- results(deseq2Data.transfected.mavs, alpha=0.05,list(c("Condition_Infected_vs_NonInfected", "TransfectionsiBrf1.ConditionInfected")), tidy=TRUE)

# View summary of results
summary(deseq2Results.transfected.mavs)

# Get GeneIDs
gns.deseq2Results.transfected.mavs <- gns
names(gns.deseq2Results.transfected.mavs) <- c("row", "symbol")
df_results.transfected.mavs <- merge(deseq2Results.transfected.mavs.tidy, gns.deseq2Results.transfected.mavs, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.transfected.mavs <- which(deseq2Results.transfected.mavs$padj < 0.001 & abs(deseq2Results.transfected.mavs$log2FoldChange) > 3) 
df.transfected.mavs <- as.data.frame(colData(deseq2Data.transfected.mavs)[,c("Condition", "Transfection")])
vst.transfected.mavs <- vst(deseq2Data.transfected.mavs)
rownames(df.transfected.mavs) <- colnames(vst.transfected.mavs)
colnames(df.transfected.mavs) <- c("Condition", "Transfection")

# get geneids 
gns.transfected.mavs <- gns %>%
  filter(Geneid %in% row.names(vst.transfected.mavs))
row.names(vst.transfected.mavs)[match(gns.transfected.mavs[,1], row.names(vst.transfected.mavs))] <- gns.transfected.mavs[,2]

pheatmap(assay(vst.transfected.mavs)[select.transfected.mavs,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.transfected.mavs )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.transfected.mavs)


## --------------------------------------------------------------------------------------------
# Filter for PKR-KO cell line
sampleData.transfected.pkr <- sampleData.transfected %>%
  filter(Cell_line == 'PKR-KO') 


# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.transfected.pkr <- sampleData.transfected.pkr$Filename
rawCounts.matrix.transfected.pkr <- as.matrix(rawCounts[,sampleIndex.transfected.pkr])
rownames(rawCounts.matrix.transfected.pkr) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.transfected.pkr) <- sampleData.transfected.pkr$Filename
keep <- c("Condition", "SampleName","Transfection")
sampleData.transfected.pkr <- sampleData.transfected.pkr[,keep]
sampleData.transfected.pkr$Condition <- factor(sampleData.transfected.pkr$Condition)
sampleData.transfected.pkr$Transfection <- factor(sampleData.transfected.pkr$Transfection)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.transfected.pkr.ordered <- rawCounts.matrix.transfected.pkr[,unique(rownames(sampleData.transfected.pkr))]
all(colnames(rawCounts.matrix.transfected.pkr.ordered) == rownames(sampleData.transfected.pkr))

# Set the reference level for each factor
sampleData.transfected.pkr$Condition <- relevel(sampleData.transfected.pkr$Condition, ref='NonInfected')
sampleData.transfected.pkr$Transfection <- relevel(sampleData.transfected.pkr$Transfection, ref='siCtrl')

#Set NA values to 0
rawCounts.matrix.transfected.pkr.ordered[is.na(rawCounts.matrix.transfected.pkr.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.transfected.pkr <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.transfected.pkr.ordered, colData=sampleData.transfected.pkr, design= ~Transfection*Condition)




## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.transfected.pkr)
plotPCA(rld, intgroup=c('Condition', 'Transfection'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.transfected.pkr)
dim(deseq2Data.transfected.pkr[rowSums(counts(deseq2Data.transfected.pkr)) > 5, ])
deseq2Data.transfected.pkr <- deseq2Data.transfected.pkr[rowSums(counts(deseq2Data.transfected.pkr)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.transfected.pkr <- DESeq(deseq2Data.transfected.pkr, parallel = TRUE)

# Extract differential expression results
deseq2Results.transfected.pkr <- results(deseq2Data.transfected.pkr, alpha=0.05,list(c("Condition_Infected_vs_NonInfected", "TransfectionsiBrf1.ConditionInfected")))
deseq2Results.transfected.pkr.tidy <- results(deseq2Data.transfected.pkr, alpha=0.05,list(c("Condition_Infected_vs_NonInfected", "TransfectionsiBrf1.ConditionInfected")), tidy=TRUE)

# View summary of results
summary(deseq2Results.transfected.pkr)

# Get GeneIDs
gns.deseq2Results.transfected.pkr <- gns
names(gns.deseq2Results.transfected.pkr) <- c("row", "symbol")
df_results.transfected.pkr <- merge(deseq2Results.transfected.pkr.tidy, gns.deseq2Results.transfected.pkr, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.transfected.pkr <- which(deseq2Results.transfected.pkr$padj < 0.05 & abs(deseq2Results.transfected.pkr$log2FoldChange) > 3) 
df.transfected.pkr <- as.data.frame(colData(deseq2Data.transfected.pkr)[,c("Condition", "Transfection")])
vst.transfected.pkr <- vst(deseq2Data.transfected.pkr)
rownames(df.transfected.pkr) <- colnames(vst.transfected.pkr)
colnames(df.transfected.pkr) <- c("Condition", "Transfection")

# get geneids
gns.transfected.pkr <- gns %>%
  filter(Geneid %in% row.names(vst.transfected.pkr))
row.names(vst.transfected.pkr)[match(gns.transfected.pkr[,1], row.names(vst.transfected.pkr))] <- gns.transfected.pkr[,2]

# Create heatmap
pheatmap(assay(vst.transfected.pkr)[select.transfected.pkr,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.transfected.pkr )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.transfected.pkr)


## --------------------------------------------------------------------------------------------
# Filter for RL-KO cell line
sampleData.transfected.rl <- sampleData.transfected %>%
  filter(Cell_line == 'RL-KO') 


# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.transfected.rl <- sampleData.transfected.rl$Filename
rawCounts.matrix.transfected.rl <- as.matrix(rawCounts[,sampleIndex.transfected.rl])
rownames(rawCounts.matrix.transfected.rl) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.transfected.rl) <- sampleData.transfected.rl$Filename
keep <- c("Condition", "SampleName","Transfection")
sampleData.transfected.rl <- sampleData.transfected.rl[,keep]
sampleData.transfected.rl$Condition <- factor(sampleData.transfected.rl$Condition)
sampleData.transfected.rl$Transfection <- factor(sampleData.transfected.rl$Transfection)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.transfected.rl.ordered <- rawCounts.matrix.transfected.rl[,unique(rownames(sampleData.transfected.rl))]
all(colnames(rawCounts.matrix.transfected.rl.ordered) == rownames(sampleData.transfected.rl))

# Set the reference level for each factor
sampleData.transfected.rl$Condition <- relevel(sampleData.transfected.rl$Condition, ref='NonInfected')
sampleData.transfected.rl$Transfection <- relevel(sampleData.transfected.rl$Transfection, ref='siCtrl')

#Set NA values to 0
rawCounts.matrix.transfected.rl.ordered[is.na(rawCounts.matrix.transfected.rl.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.transfected.rl <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.transfected.rl.ordered, colData=sampleData.transfected.rl, design= ~Transfection*Condition)



## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.transfected.rl)
plotPCA(rld, intgroup=c('Condition', 'Transfection'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.transfected.rl)
dim(deseq2Data.transfected.rl[rowSums(counts(deseq2Data.transfected.rl)) > 5, ])
deseq2Data.transfected.rl <- deseq2Data.transfected.rl[rowSums(counts(deseq2Data.transfected.rl)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.transfected.rl <- DESeq(deseq2Data.transfected.rl, parallel = TRUE)

# Extract differential expression results
deseq2Results.transfected.rl <- results(deseq2Data.transfected.rl, alpha=0.05,list(c("Condition_Infected_vs_NonInfected", "TransfectionsiBrf1.ConditionInfected")))
deseq2Results.transfected.rl.tidy  <- results(deseq2Data.transfected.rl, alpha=0.05,list(c("Condition_Infected_vs_NonInfected", "TransfectionsiBrf1.ConditionInfected")), tidy=TRUE)

# View summary of results
summary(deseq2Results.transfected.rl)

# Get GeneIDs
gns.deseq2Results.transfected.rl <- gns
names(gns.deseq2Results.transfected.rl) <- c("row", "symbol")
df_results.transfected.rl <- merge(deseq2Results.transfected.rl.tidy, gns.deseq2Results.transfected.rl, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.transfected.rl <- which(deseq2Results.transfected.rl$padj < 0.05 & abs(deseq2Results.transfected.rl$log2FoldChange) > 3.3)
df.transfected.rl <- as.data.frame(colData(deseq2Data.transfected.rl)[,c("Condition", "Transfection")])
vst.transfected.rl <- vst(deseq2Data.transfected.rl)
rownames(df.transfected.rl) <- colnames(vst.transfected.rl)
colnames(df.transfected.rl) <- c("Condition", "Transfection")

# get geneids
gns.transfected.rl <- gns %>%
  filter(Geneid %in% row.names(vst.transfected.rl))
row.names(vst.transfected.rl)[match(gns.transfected.rl[,1], row.names(vst.transfected.rl))] <- gns.transfected.rl[,2]

# Create Heatmap
pheatmap(assay(vst.transfected.rl)[select.transfected.rl,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.transfected.rl )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.transfected.rl)


## --------------------------------------------------------------------------------------------
# Filter for TLR3-KO cell line
sampleData.transfected.tlr3 <- sampleData.transfected %>%
  filter(Cell_line == 'TLR3-KO') 


# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.transfected.tlr3 <- sampleData.transfected.tlr3$Filename
rawCounts.matrix.transfected.tlr3 <- as.matrix(rawCounts[,sampleIndex.transfected.tlr3])
rownames(rawCounts.matrix.transfected.tlr3) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.transfected.tlr3) <- sampleData.transfected.tlr3$Filename
keep <- c("Condition", "SampleName","Transfection")
sampleData.transfected.tlr3 <- sampleData.transfected.tlr3[,keep]
sampleData.transfected.tlr3$Condition <- factor(sampleData.transfected.tlr3$Condition)
sampleData.transfected.tlr3$Transfection <- factor(sampleData.transfected.tlr3$Transfection)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.transfected.tlr3.ordered <- rawCounts.matrix.transfected.tlr3[,unique(rownames(sampleData.transfected.tlr3))]
all(colnames(rawCounts.matrix.transfected.tlr3.ordered) == rownames(sampleData.transfected.tlr3))

# Set the reference level for each factor
sampleData.transfected.tlr3$Condition <- relevel(sampleData.transfected.tlr3$Condition, ref='NonInfected')
sampleData.transfected.tlr3$Transfection <- relevel(sampleData.transfected.tlr3$Transfection, ref='siCtrl')

#Set NA values to 0
rawCounts.matrix.transfected.tlr3.ordered[is.na(rawCounts.matrix.transfected.tlr3.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.transfected.tlr3 <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.transfected.tlr3.ordered, colData=sampleData.transfected.tlr3, design= ~Transfection*Condition)




## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.transfected.tlr3)
plotPCA(rld, intgroup=c('Condition', 'Transfection'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.transfected.tlr3)
dim(deseq2Data.transfected.tlr3[rowSums(counts(deseq2Data.transfected.tlr3)) > 5, ])
deseq2Data.transfected.tlr3 <- deseq2Data.transfected.tlr3[rowSums(counts(deseq2Data.transfected.tlr3)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.transfected.tlr3 <- DESeq(deseq2Data.transfected.tlr3, parallel = TRUE)

# Extract differential expression results
deseq2Results.transfected.tlr3 <- results(deseq2Data.transfected.tlr3, alpha=0.05,list(c("Condition_Infected_vs_NonInfected", "TransfectionsiBrf1.ConditionInfected")))
deseq2Results.transfected.tlr3.tidy  <- results(deseq2Data.transfected.tlr3, alpha=0.05,list(c("Condition_Infected_vs_NonInfected", "TransfectionsiBrf1.ConditionInfected")), tidy=TRUE)

# View summary of results
summary(deseq2Results.transfected.tlr3)

# Get GeneIDs
gns.deseq2Results.transfected.tlr3 <- gns
names(gns.deseq2Results.transfected.tlr3) <- c("row", "symbol")
df_results.transfected.tlr3 <- merge(deseq2Results.transfected.tlr3.tidy, gns.deseq2Results.transfected.tlr3, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.transfected.tlr3 <- which(deseq2Results.transfected.tlr3$padj < 0.0005 & abs(deseq2Results.transfected.tlr3$log2FoldChange) > 3.2)
df.transfected.tlr3 <- as.data.frame(colData(deseq2Data.transfected.tlr3)[,c("Condition", "Transfection")])
vst.transfected.tlr3 <- vst(deseq2Data.transfected.tlr3)
rownames(df.transfected.tlr3) <- colnames(vst.transfected.tlr3)
colnames(df.transfected.tlr3) <- c("Condition", "Transfection")

# get geneids
gns.transfected.tlr3 <- gns %>%
  filter(Geneid %in% row.names(vst.transfected.tlr3))
row.names(vst.transfected.tlr3)[match(gns.transfected.tlr3[,1], row.names(vst.transfected.tlr3))] <- gns.transfected.tlr3[,2]

# Create Heatmap
pheatmap(assay(vst.transfected.tlr3)[select.transfected.tlr3,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.transfected.tlr3 )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.transfected.tlr3)


## --------------------------------------------------------------------------------------------
#Filter for only infected cells
sampleData.infected <- sampleData %>%
  filter(Condition != 'NonInfected')


## --------------------------------------------------------------------------------------------
# Filter for WT cell line
sampleData.infected.wt <- sampleData.infected %>%
  filter(Cell_line == 'A549') %>%
  filter(Batch == 'nR052')

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.infected.wt <- sampleData.infected.wt$Filename
rawCounts.matrix.infected.wt <- as.matrix(rawCounts[,sampleIndex.infected.wt])
rownames(rawCounts.matrix.infected.wt) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.infected.wt) <- sampleData.infected.wt$Filename
keep <- c("Condition", "SampleName","Transfection", "Binary_Transfection")
sampleData.infected.wt <- sampleData.infected.wt[,keep]
sampleData.infected.wt$Binary_Transfection <- factor(sampleData.infected.wt$Binary_Transfection)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.infected.wt.ordered <- rawCounts.matrix.infected.wt[,unique(rownames(sampleData.infected.wt))]
all(colnames(rawCounts.matrix.infected.wt.ordered) == rownames(sampleData.infected.wt))

# Set the reference level for each factor
sampleData.infected.wt$Binary_Transfection <- relevel(sampleData.infected.wt$Binary_Transfection, ref='untransfected')

#Set NA values to 0
rawCounts.matrix.infected.wt.ordered[is.na(rawCounts.matrix.infected.wt.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.infected.wt <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.infected.wt.ordered, colData=sampleData.infected.wt, design= ~Binary_Transfection)


## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.infected.wt)
plotPCA(rld, intgroup=c('Binary_Transfection'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.infected.wt)
dim(deseq2Data.infected.wt[rowSums(counts(deseq2Data.infected.wt)) > 5, ])
deseq2Data.infected.wt <- deseq2Data.infected.wt[rowSums(counts(deseq2Data.infected.wt)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.infected.wt <- DESeq(deseq2Data.infected.wt, parallel = TRUE)

# Extract differential expression results
deseq2Results.infected.wt <- results(deseq2Data.infected.wt, alpha=0.05, contrast=c("Binary_Transfection", "transfected", "untransfected"))
deseq2Results.infected.wt.tidy <- results(deseq2Data.infected.wt, alpha=0.05, contrast=c("Binary_Transfection", "transfected", "untransfected"), tidy=TRUE)

# View summary of results
summary(deseq2Results.infected.wt)

# GeneIDs retrieved February 2022
gns.deseq2Results.infected.wt <- gns
names(gns.deseq2Results.infected.wt) <- c("row", "symbol")
df_results.infected.wt <- merge(deseq2Results.infected.wt.tidy, gns.deseq2Results.infected.wt, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.infected.wt <- which(deseq2Results.infected.wt$padj < 0.000000005)
df.infected.wt <- as.data.frame(colData(deseq2Data.infected.wt)[,c("Binary_Transfection")])
vst.infected.wt <- vst(deseq2Data.infected.wt)
rownames(df.infected.wt) <- colnames(vst.infected.wt)
colnames(df.infected.wt) <- c("Binary_Transfection")

# get geneids
gns.infected.wt <- gns %>%
  filter(Geneid %in% row.names(vst.infected.wt))
row.names(vst.infected.wt)[match(gns.infected.wt[,1], row.names(vst.infected.wt))] <- gns.infected.wt[,2]

# Create Heatmap
pheatmap(assay(vst.infected.wt)[select.infected.wt,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.infected.wt )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.infected.wt)


## --------------------------------------------------------------------------------------------
# Filter for WT cell line
sampleData.infected.mavs <- sampleData.infected %>%
  filter(Cell_line == 'MAVS-KO') 

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.infected.mavs <- sampleData.infected.mavs$Filename
rawCounts.matrix.infected.mavs <- as.matrix(rawCounts[,sampleIndex.infected.mavs])
rownames(rawCounts.matrix.infected.mavs) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.infected.mavs) <- sampleData.infected.mavs$Filename
keep <- c("Condition", "SampleName","Transfection", "Binary_Transfection")
sampleData.infected.mavs <- sampleData.infected.mavs[,keep]
sampleData.infected.mavs$Condition <- factor(sampleData.infected.mavs$Condition)
sampleData.infected.mavs$Transfection <- factor(sampleData.infected.mavs$Transfection)
sampleData.infected.mavs$Binary_Transfection <- factor(sampleData.infected.mavs$Binary_Transfection)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.infected.mavs.ordered <- rawCounts.matrix.infected.mavs[,unique(rownames(sampleData.infected.mavs))]
all(colnames(rawCounts.matrix.infected.mavs.ordered) == rownames(sampleData.infected.mavs))

# Set the reference level for each factor
sampleData.infected.mavs$Transfection <- relevel(sampleData.infected.mavs$Transfection, ref='siCtrl')
sampleData.infected.mavs$Binary_Transfection <- relevel(sampleData.infected.mavs$Binary_Transfection, ref='untransfected')

#Set NA values to 0
rawCounts.matrix.infected.mavs.ordered[is.na(rawCounts.matrix.infected.mavs.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.infected.mavs <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.infected.mavs.ordered, colData=sampleData.infected.mavs, design= ~Binary_Transfection)




## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.infected.mavs)
plotPCA(rld, intgroup=c('Binary_Transfection'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.infected.mavs)
dim(deseq2Data.infected.mavs[rowSums(counts(deseq2Data.infected.mavs)) > 5, ])
deseq2Data.infected.mavs <- deseq2Data.infected.mavs[rowSums(counts(deseq2Data.infected.mavs)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.infected.mavs <- DESeq(deseq2Data.infected.mavs, parallel = FALSE)

# Extract differential expression results
deseq2Results.infected.mavs <- results(deseq2Data.infected.mavs, alpha=0.05, contrast=c("Binary_Transfection", "transfected", "untransfected"))
deseq2Results.infected.mavs.tidy <- results(deseq2Data.infected.mavs, alpha=0.05, contrast=c("Binary_Transfection", "transfected", "untransfected"), tidy=TRUE)

# View summary of results
summary(deseq2Results.infected.mavs)

# Get GeneIDs
gns.deseq2Results.infected.mavs <- gns
names(gns.deseq2Results.infected.mavs) <- c("row", "symbol")
df_results.infected.mavs <- merge(deseq2Results.infected.mavs.tidy, gns.deseq2Results.infected.mavs, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.infected.mavs <- which(deseq2Results.infected.mavs$padj < 0.5)
df.infected.mavs <- as.data.frame(colData(deseq2Data.infected.mavs)[,c("Binary_Transfection")])
vst.infected.mavs <- vst(deseq2Data.infected.mavs)
rownames(df.infected.mavs) <- colnames(vst.infected.mavs)
colnames(df.infected.mavs) <- c("Binary_Transfection")

# get geneids 
gns.infected.mavs <- gns %>%
  filter(Geneid %in% row.names(vst.infected.mavs))
row.names(vst.infected.mavs)[match(gns.infected.mavs[,1], row.names(vst.infected.mavs))] <- gns.infected.mavs[,2]

# Create heatmap
pheatmap(assay(vst.infected.mavs)[select.infected.mavs,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.infected.mavs )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.infected.mavs)


## --------------------------------------------------------------------------------------------
# Filter for PKR cell line
sampleData.infected.pkr <- sampleData.infected %>%
  filter(Cell_line == 'PKR-KO')

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex <- sampleData.infected.pkr$Filename
rawCounts.matrix.infected.pkr <- as.matrix(rawCounts[,sampleIndex])
rownames(rawCounts.matrix.infected.pkr) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.infected.pkr) <- sampleData.infected.pkr$Filename
keep <- c("Condition", "SampleName","Transfection", "Binary_Transfection")
sampleData.infected.pkr <- sampleData.infected.pkr[,keep]
sampleData.infected.pkr$Condition <- factor(sampleData.infected.pkr$Condition)
sampleData.infected.pkr$Transfection <- factor(sampleData.infected.pkr$Transfection)
sampleData.infected.pkr$Binary_Transfection <- factor(sampleData.infected.pkr$Binary_Transfection)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.infected.pkr.ordered <- rawCounts.matrix.infected.pkr[,unique(rownames(sampleData.infected.pkr))]
all(colnames(rawCounts.matrix.infected.pkr.ordered) == rownames(sampleData.infected.pkr))

# Set the reference level for each factor
sampleData.infected.pkr$Transfection <- relevel(sampleData.infected.pkr$Transfection, ref='siCtrl')
sampleData.infected.pkr$Binary_Transfection <- relevel(sampleData.infected.pkr$Binary_Transfection, ref='untransfected')

#Set NA values to 0
rawCounts.matrix.infected.pkr.ordered[is.na(rawCounts.matrix.infected.pkr.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.infected.pkr <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.infected.pkr.ordered, colData=sampleData.infected.pkr, design= ~Binary_Transfection)




## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.infected.pkr)
plotPCA(rld, intgroup=c('Binary_Transfection'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.infected.pkr)
dim(deseq2Data.infected.pkr[rowSums(counts(deseq2Data.infected.pkr)) > 5, ])
deseq2Data.infected.pkr <- deseq2Data.infected.pkr[rowSums(counts(deseq2Data.infected.pkr)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.infected.pkr <- DESeq(deseq2Data.infected.pkr, parallel = FALSE)

# Extract differential expression results
deseq2Results.infected.pkr <- results(deseq2Data.infected.pkr, alpha=0.05, contrast=c("Binary_Transfection", "transfected", "untransfected"))
deseq2Results.infected.pkr.tidy <- results(deseq2Data.infected.pkr, alpha=0.05, contrast=c("Binary_Transfection", "transfected", "untransfected"), tidy=TRUE)

# View summary of results
summary(deseq2Results.infected.pkr)

# Get GeneIDs
gns.deseq2Results.infected.pkr <- gns
names(gns.deseq2Results.infected.pkr) <- c("row", "symbol")
df_results.infected.pkr <- merge(deseq2Results.infected.pkr.tidy, gns.deseq2Results.infected.pkr, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.infected.pkr <- which(deseq2Results.infected.pkr$padj < 0.05)
df.infected.pkr <- as.data.frame(colData(deseq2Data.infected.pkr)[,c("Binary_Transfection")])
vst.infected.pkr <- vst(deseq2Data.infected.pkr)
rownames(df.infected.pkr) <- colnames(vst.infected.pkr)
colnames(df.infected.pkr) <- c("Binary_Transfection")

# get geneids
gns.infected.pkr <- gns %>%
  filter(Geneid %in% row.names(vst.infected.pkr))
row.names(vst.infected.pkr)[match(gns.infected.pkr[,1], row.names(vst.infected.pkr))] <- gns.infected.pkr[,2]

# Create Heatmap
pheatmap(assay(vst.infected.pkr)[select.infected.pkr,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.infected.pkr )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.infected.pkr)


## --------------------------------------------------------------------------------------------
# Filter for WT cell line
sampleData.infected.rl <- sampleData.infected %>%
  filter(Cell_line == 'RL-KO')

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.infected.rl <- sampleData.infected.rl$Filename
rawCounts.matrix.infected.rl <- as.matrix(rawCounts[,sampleIndex.infected.rl])
rownames(rawCounts.matrix.infected.rl) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.infected.rl) <- sampleData.infected.rl$Filename
keep <- c("Condition", "SampleName","Transfection", "Binary_Transfection")
sampleData.infected.rl <- sampleData.infected.rl[,keep]
sampleData.infected.rl$Condition <- factor(sampleData.infected.rl$Condition)
sampleData.infected.rl$Transfection <- factor(sampleData.infected.rl$Transfection)
sampleData.infected.rl$Binary_Transfection <- factor(sampleData.infected.rl$Binary_Transfection)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.infected.rl.ordered <- rawCounts.matrix.infected.rl[,unique(rownames(sampleData.infected.rl))]
all(colnames(rawCounts.matrix.infected.rl.ordered) == rownames(sampleData.infected.rl))

# Set the reference level for each factor
sampleData.infected.rl$Transfection <- relevel(sampleData.infected.rl$Transfection, ref='siCtrl')
sampleData.infected.rl$Binary_Transfection <- relevel(sampleData.infected.rl$Binary_Transfection, ref='untransfected')

#Set NA values to 0
rawCounts.matrix.infected.rl.ordered[is.na(rawCounts.matrix.infected.rl.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.infected.rl <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.infected.rl.ordered, colData=sampleData.infected.rl, design= ~Binary_Transfection)




## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.infected.rl)
plotPCA(rld, intgroup=c('Binary_Transfection'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.infected.rl)
dim(deseq2Data.infected.rl[rowSums(counts(deseq2Data.infected.rl)) > 5, ])
deseq2Data.infected.rl <- deseq2Data.infected.rl[rowSums(counts(deseq2Data.infected.rl)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.infected.rl <- DESeq(deseq2Data.infected.rl, parallel = FALSE)

# Extract differential expression results
deseq2Results.infected.rl <- results(deseq2Data.infected.rl, alpha=0.05, contrast=c("Binary_Transfection", "transfected", "untransfected"))
deseq2Results.infected.rl.tidy <- results(deseq2Data.infected.rl, alpha=0.05, contrast=c("Binary_Transfection", "transfected", "untransfected"), tidy=TRUE)

# View summary of results
summary(deseq2Results.infected.rl)

# get geneIDs
gns.deseq2Results.infected.rl <- gns
names(gns.deseq2Results.infected.rl) <- c("row", "symbol")
df_results.infected.rl <- merge(deseq2Results.infected.rl.tidy, gns.deseq2Results.infected.rl, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.infected.rl <- which(deseq2Results.infected.rl$padj < 0.0000000005)
df.infected.rl <- as.data.frame(colData(deseq2Data.infected.rl)[,c("Binary_Transfection")])
vst.infected.rl <- vst(deseq2Data.infected.rl)
rownames(df.infected.rl) <- colnames(vst.infected.rl)
colnames(df.infected.rl) <- c("Binary_Transfection")

# get geneids
gns.infected.rl <- gns %>%
  filter(Geneid %in% row.names(vst.infected.rl))
row.names(vst.infected.rl)[match(gns.infected.rl[,1], row.names(vst.infected.rl))] <- gns.infected.rl[,2]

# Create Heatmap
pheatmap(assay(vst.infected.rl)[select.infected.rl,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.infected.rl )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.infected.rl)


## --------------------------------------------------------------------------------------------
# Filter for tlr3 cell line
sampleData.infected.tlr3 <- sampleData.infected %>%
  filter(Cell_line == 'TLR3-KO')

# Convert count data to a matrix of appropriate form that DEseq2 can read
geneID <- rawCounts$Geneid
sampleIndex.infected.tlr3 <- sampleData.infected.tlr3$Filename
rawCounts.matrix.infected.tlr3 <- as.matrix(rawCounts[,sampleIndex.infected.tlr3])
rownames(rawCounts.matrix.infected.tlr3) <- geneID

# Convert sample variable mappings to an appropriate form that DESeq2 can read
rownames(sampleData.infected.tlr3) <- sampleData.infected.tlr3$Filename
keep <- c("Condition", "SampleName","Transfection", "Binary_Transfection")
sampleData.infected.tlr3 <- sampleData.infected.tlr3[,keep]
sampleData.infected.tlr3$Condition <- factor(sampleData.infected.tlr3$Condition)
sampleData.infected.tlr3$Transfection <- factor(sampleData.infected.tlr3$Transfection)
sampleData.infected.tlr3$Binary_Transfection <- factor(sampleData.infected.tlr3$Binary_Transfection)

# Put the columns of the count data in the same order as rows names of the sample mapping, then make sure it worked
rawCounts.matrix.infected.tlr3.ordered <- rawCounts.matrix.infected.tlr3[,unique(rownames(sampleData.infected.tlr3))]
all(colnames(rawCounts.matrix.infected.tlr3.ordered) == rownames(sampleData.infected.tlr3))

# Set the reference level for each factor
sampleData.infected.tlr3$Transfection <- relevel(sampleData.infected.tlr3$Transfection, ref='siCtrl')
sampleData.infected.tlr3$Binary_Transfection <- relevel(sampleData.infected.tlr3$Binary_Transfection, ref='untransfected')

#Set NA values to 0
rawCounts.matrix.infected.tlr3.ordered[is.na(rawCounts.matrix.infected.tlr3.ordered)] <- 0

# Create the DEseq2DataSet object
deseq2Data.infected.tlr3 <- DESeqDataSetFromMatrix(countData=rawCounts.matrix.infected.tlr3.ordered, colData=sampleData.infected.tlr3, design= ~Binary_Transfection)


## --------------------------------------------------------------------------------------------
# Plot PCA of deseq2Data
rld <- rlog(deseq2Data.infected.tlr3)
plotPCA(rld, intgroup=c('Binary_Transfection'))


## --------------------------------------------------------------------------------------------
# Perform pre-filtering of the data
dim(deseq2Data.infected.tlr3)
dim(deseq2Data.infected.tlr3[rowSums(counts(deseq2Data.infected.tlr3)) > 5, ])
deseq2Data.infected.tlr3 <- deseq2Data.infected.tlr3[rowSums(counts(deseq2Data.infected.tlr3)) > 5, ]


# Register the number of cores to use
register(MulticoreParam(4))

# Run pipeline for differential expression steps (if you set up parallel processing, set parallel = TRUE here)
deseq2Data.infected.tlr3 <- DESeq(deseq2Data.infected.tlr3, parallel = FALSE)

# Extract differential expression results
deseq2Results.infected.tlr3 <- results(deseq2Data.infected.tlr3, alpha=0.05, contrast=c("Binary_Transfection", "transfected", "untransfected"))
deseq2Results.infected.tlr3.tidy <- results(deseq2Data.infected.tlr3, alpha=0.05, contrast=c("Binary_Transfection", "transfected", "untransfected"), tidy=TRUE)

# View summary of results
summary(deseq2Results.infected.tlr3)

#Get GeneIDs
gns.deseq2Results.infected.tlr3 <- gns
names(gns.deseq2Results.infected.tlr3) <- c("row", "symbol")
df_results.infected.tlr3 <- merge(deseq2Results.infected.tlr3.tidy, gns.deseq2Results.infected.tlr3, all.x=TRUE)



## --------------------------------------------------------------------------------------------
# Extract top genes with significant padj
select.infected.tlr3 <- which(deseq2Results.infected.tlr3$padj < 0.005)
df.infected.tlr3 <- as.data.frame(colData(deseq2Data.infected.tlr3)[,c("Binary_Transfection")])
vst.infected.tlr3 <- vst(deseq2Data.infected.tlr3)
rownames(df.infected.tlr3) <- colnames(vst.infected.tlr3)
colnames(df.infected.tlr3) <- c("Binary_Transfection")

# get geneids
gns.infected.tlr3 <- gns %>%
  filter(Geneid %in% row.names(vst.infected.tlr3))
row.names(vst.infected.tlr3)[match(gns.infected.tlr3[,1], row.names(vst.infected.tlr3))] <- gns.infected.tlr3[,2]

# Create Heatmap
pheatmap(assay(vst.infected.tlr3)[select.infected.tlr3,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col =df.infected.tlr3 )


## --------------------------------------------------------------------------------------------
# Using DEseq2 built in method
plotMA(deseq2Results.infected.tlr3)

