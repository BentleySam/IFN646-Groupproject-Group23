library(readr)
library (DESeq2)
library (tidyverse)

counts <- read.table("GSE159717_rnaseq_deseq_5dpi_counts_raw.tsv",
                     header = TRUE, row.names = 1, sep = "\t")

head(counts)
dim(counts)

counts_filtered <- counts[counts$gene_biotype == "protein_coding", ]

dim (counts_filtered)

counts_filtered <- counts_filtered[, c("S_2_Rem_5dpi_S70001",  "S_2_SARS_5dpi_S70003", "S_2_mock_5dpi_S70002",
                                       "S_3_Rem_5dpi_S69995",  "S_3_SARS_5dpi_S69996", "S_3_mock_5dpi_S69997")]

sample_names <- colnames(counts_filtered)

coldata <- data.frame(
  condition = dplyr::case_when(
    grepl("mock", sample_names) ~ "mock",
    grepl("SARS", sample_names) ~ "SARS",
    grepl("Rem",  sample_names) ~ "Rem"
  ),
  row.names = sample_names
)

all(colnames(counts_filtered) %in% rownames(coldata))

all(colnames(counts_filtered) == rownames(coldata))

dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                       colData = coldata,
                       design = ~ condition)

dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]

dds <- DESeq(dds)

res = results(dds, contrast=c("condition", "mock", "SARS"), alpha = 1e-5)
res

res_ordered <- res[order(res$padj), ]
res_ordered

sig_res <- subset(res_ordered, padj < 0.05 & abs(log2FoldChange) > 1.0)
sig_res

plotMA(res)

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(res, lab=rownames(res), x = 'log2FoldChange', y='pvalue')


vsd <- vst (dds)

plotPCA(vsd, intgroup = c("condition"))

pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentvar <- round(100 * attr(pcaData, "% Var"))

rld <- rlog(dds, blind=TRUE)

plotPCA(rld, intgroup="condition")

rld <- rlog(dds, blind=T)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))

df <- cbind(coldata, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = condition))


#--
#EdgeR  
library(edgeR)

counts <- read.table("GSE159717_rnaseq_deseq_5dpi_counts_raw.tsv",
                     header = TRUE, row.names = 1, sep = "\t")

# Filter protein-coding genes
counts_filtered2 <- counts[counts$gene_biotype == "protein_coding", ]

# Select columns you want
counts_filtered2 <- counts_filtered2[, c("S_2_Rem_5dpi_S70001",  "S_2_SARS_5dpi_S70003", "S_2_mock_5dpi_S70002",
                                       "S_3_Rem_5dpi_S69995",  "S_3_SARS_5dpi_S69996", "S_3_mock_5dpi_S69997")]

sample_names <- colnames(counts_filtered2)

group <- factor(dplyr::case_when(
  grepl("mock", sample_names, ignore.case = TRUE) ~ "mock",
  grepl("SARS", sample_names, ignore.case = TRUE) ~ "SARS",
  grepl("Rem",  sample_names, ignore.case = TRUE) ~ "Rem"
))

coldata2 <- data.frame(condition = factor(group), row.names = sample_names)


dge <- DGEList(counts = counts_filtered2, group = coldata2$condition)

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

dge <- estimateDisp(dge)

fit <- glmFit(dge, design = model.matrix(~condition, data=coldata2))
lrt <- glmLRT(fit, contrast = c(0,1,0))

topTags(lrt)

res2 <- as.data.frame(lrt$table)
res2$FDR <- p.adjust(res2$PValue, method="BH")

sig_res2 <- subset(res2, FDR < 0.1 & abs(logFC) > 0.5)
sig_res2

# Volcano plot
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'logFC',
                y = 'PValue')

logCPM <- cpm(dge, log=TRUE)
pca <- prcomp(t(logCPM))

pca_df <- cbind(coldata2, pca$x)
ggplot(pca_df, aes(x=PC1, y=PC2, color=condition)) + geom_point(size=3)


#-----
library(limma)

counts <- read.table("GSE159717_rnaseq_deseq_5dpi_counts_raw.tsv",
                     header = TRUE, row.names = 1, sep = "\t")

counts_filtered3 <- counts[counts$gene_biotype == "protein_coding", ]

# Select samples
counts_filtered3 <- counts_filtered3[, c("S_2_Rem_5dpi_S70001",  "S_2_SARS_5dpi_S70003", "S_2_mock_5dpi_S70002",
                                       "S_3_Rem_5dpi_S69995",  "S_3_SARS_5dpi_S69996", "S_3_mock_5dpi_S69997")]

# Define group/condition
sample_names <- colnames(counts_filtered3)
group <- factor(dplyr::case_when(
  grepl("mock", sample_names, ignore.case = TRUE) ~ "mock",
  grepl("SARS", sample_names, ignore.case = TRUE) ~ "SARS",
  grepl("Rem",  sample_names, ignore.case = TRUE) ~ "Rem"
))
coldata3 <- data.frame(condition = group, row.names = sample_names)

dge2 <- DGEList(counts = counts_filtered3, group = coldata3$condition)
keep <- filterByExpr(dge2)
dge2 <- dge2[keep, , keep.lib.sizes = FALSE]

dge2 <- calcNormFactors(dge2)

design <- model.matrix(~condition, data=coldata3)

v <- voom(dge2, design, plot=TRUE)  

fit <- lmFit(v, design)
fit <- eBayes(fit)  

res3 <- topTable(fit, coef="conditionSARS", number=Inf, sort.by="P")

res3$FDR <- p.adjust(res3$P.Value, method="BH")

sig_res3 <- subset(res3, FDR < 0.05 & abs(logFC) > 1)
sig_res3

EnhancedVolcano(res3,
                lab = rownames(res3),
                x = 'logFC',
                y = 'P.Value')

