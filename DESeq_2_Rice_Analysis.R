library(stringr)
library(DESeq2)
library(ggplot2)
library(ashr)
library(tidyverse)
library("dplyr")
library(apeglm)
library(vsn)
library(gplots)

getwd()
list.files()
setwd("E:/projects/rice/pipeline_result/heat_drought/data/workflow_SE/results/featureCounts")
counts <- read.table("counts.txt", sep = "", head = T, skip = 1, row.names = "Geneid")
countsnew <- counts[(-c(1:5))]
colnames(countsnew)= str_split_fixed(colnames(countsnew), "\\.",6)[,1]

mycols <- data.frame(row.names = (colnames(countsnew)))
coldata <- data.frame(mycols, condition = factor(c(rep("combined", 3), rep("control", 3))))
dds = DESeqDataSetFromMatrix(countData = countsnew, colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10,]
dds <- DESeq(dds)

vst <- vst(dds, blind = FALSE)

plotPCA(vst, intgroup="condition", ntop=nrow(counts(dds)))

a <- DESeq2::plotPCA(vst, intgroup="condition")
a = geom_label(aes(label = coldata$condition),)
nudge <- position_nudge(y = 1)
boxplot(assay(vst), col = c("Red", "Red", "Red", "Green", "Green", "Green"), pch = ".",
        vertical = TRUE, cex.axis=0.5, main = "Boxplot of heat and drought affected rice using vst method",
        las = 2, ylab = "assay(vst)", xlab = "Samples", ylim=c(-10, 30),
        font.main = 5, font.axis = 0.5, font.lab=2 )
cU <- cor(as.matrix(assay(vst)))
cols <- c("dodgerblue3", "firebrick3")[coldata$condition]
heatmap.2(cU, symm = TRUE, col = colorRampPalette(c("darkblue", "white"))(100),
          labCol = colnames(cU), labRow = colnames(cU),
          distfun = function(c) as.dist(1 - c),
          trace = "none",
          Colv = TRUE, cexRow = 0.9, cexCol = 0.9, key = F,
          font=2,
          RowSideColors = cols, ColSideColors = cols)
plotDispEsts(dds)
res <- results(dds, contrast = c("condition", "combined", "control"))
summary(res)
grp.mean <- sapply(levels(dds$condition),
                   function(lvl)
                     rowMeans(counts(dds, normalized=TRUE)[,dds$condition== lvl]))
norm.counts <- counts(dds, normalized=TRUE)
all <- data.frame(res, assay(vst))
nrow(all)
write.table(all, file = "rice_combined_stress_main.csv", sep = ",")
sum(res05$padj < 0.01, na.rm=TRUE)
table(res05$padj < 0.1 & res05$log2FoldChange > 0)
table(res05$padj < 0.01 & res05$log2FoldChange < 0)
table(res05$padj < 0.01 & res05$log2FoldChange < 0)