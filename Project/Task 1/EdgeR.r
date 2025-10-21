library(edgeR); library(limma); library(readr); library(dplyr); library(tibble)


cts <- as.matrix(read_tsv("counts.tsv"))
smp <- read_csv("samples.csv")
stopifnot(all(colnames(cts) == smp$sample))
smp$condition <- factor(smp$condition, levels = c("Control","COVID"))


# Filter: keep genes with counts >=10 in ≥20% of samples
keep <- rowSums(cts >= 10) >= max(2, floor(ncol(cts)*0.2))
cts  <- cts[keep,]


# Design (add covariates if present, e.g. ~ batch + age + sex + condition)
design <- model.matrix(~ smp$condition)

# ---- edgeR ----
dge <- DGEList(cts); dge <- calcNormFactors(dge, "TMM")
dge <- estimateDisp(dge, design)
fit  <- glmFit(dge, design)
lrt  <- glmLRT(fit, coef = "smp$conditionCOVID")
res_edgeR <- topTags(lrt, n=Inf)$table |>
  rownames_to_column("gene") |>
  transmute(gene, logFC_edgeR = logFC, p_edgeR = PValue, FDR_edgeR = FDR)

# ---- limma-voom ----
v <- voom(dge, design, plot=FALSE)
fit2 <- eBayes(lmFit(v, design))
res_limma <- topTable(fit2, coef = "smp$conditionCOVID", number = Inf, sort.by="none") |>
  rownames_to_column("gene") |>
  transmute(gene, logFC_limma = logFC, p_limma = P.Value, FDR_limma = adj.P.Val)


res <- inner_join(res_edgeR, res_limma, by="gene") |>
  mutate(
    agree_dir = sign(logFC_edgeR) == sign(logFC_limma),
    pass_edgeR = FDR_edgeR < 0.05 & abs(logFC_edgeR) >= 1,
    pass_limma = FDR_limma < 0.05 & abs(logFC_limma) >= 1,
    consensus_strict = pass_edgeR & pass_limma & agree_dir
  )
write_csv(res, "DE_consensus_edgeR_limma.csv")
write_csv(filter(res, consensus_strict), "DE_consensus_STRICT.csv")