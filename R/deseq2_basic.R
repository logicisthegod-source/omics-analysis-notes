# Basic DESeq2 workflow for RNA-seq differential expression
set.seed(42)

library(DESeq2)
library(tidyverse)

# Load count matrix and sample metadata
counts <- read.csv("data/counts.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("data/samples.csv", row.names = 1)

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = coldata,
  design = ~ condition
)

# Pre-filter low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run differential expression analysis
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treatment", "control"))

# Shrink log2 fold changes
resLFC <- lfcShrink(dds, coef = 2, res = res, type = "apeglm")

# Significant genes: padj < 0.05 and |log2FC| > 1
sig <- subset(as.data.frame(resLFC), padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig, "results/deseq2_significant.csv")

summary(res)
