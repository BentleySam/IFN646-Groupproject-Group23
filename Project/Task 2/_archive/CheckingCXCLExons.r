library(EnsDb.Hsapiens.v86)
library(AnnotationFilter)
library(GenomeInfoDb)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
Hs  <- BSgenome.Hsapiens.UCSC.hg38
edb <- EnsDb.Hsapiens.v86
gene_symbol <- "CXCL10"

gr_gene <- genes(edb, filter = SymbolFilter(gene_symbol))
ensg_id <- gr_gene$gene_id[1]
txs <- transcripts(edb, filter = GeneIdFilter(ensg_id))
cat("N transcripts:", length(txs), "\n")

ex_tx <- exons(edb, filter = TxIdFilter(txs$tx_id), columns = c("tx_id","exon_idx"))
first_exons_tx <- ex_tx[ex_tx$exon_idx == 1]
first_exons_tx <- keepStandardChromosomes(first_exons_tx, pruning.mode = "coarse")
seqlevelsStyle(first_exons_tx) <- seqlevelsStyle(Hs)

df <- data.frame(
  tx_id = first_exons_tx$tx_id,
  chr   = as.character(seqnames(first_exons_tx)),
  start = start(first_exons_tx),
  end   = end(first_exons_tx),
  strand= as.character(strand(first_exons_tx))
)
cat("Unique first-exon loci:\n")
print(unique(df[ ,c("chr","start","end","strand")]))


library(Biostrings)

ex_all <- exonsBy(edb, by = "gene")[[ensg_id]]
ex_all <- keepStandardChromosomes(ex_all, pruning.mode = "coarse")
seqlevelsStyle(ex_all) <- seqlevelsStyle(Hs)

seqs <- getSeq(Hs, ex_all)
names(seqs) <- paste0(
  gene_symbol, "|exon:", seq_along(ex_all), "|",
  as.character(seqnames(ex_all)), ":", start(ex_all), "-", end(ex_all),
  "(", as.character(strand(ex_all)), ")"
)

out_fa <- "~/crackling_inputs/CXCL10_all_exonsv2.fa"
writeXStringSet(seqs, out_fa)
cat("Wrote", length(seqs), "exon sequences to", out_fa, "\n")