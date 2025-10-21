suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
})



# If conflicted is active, prefer base::intersect to avoid errors later
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(base::intersect)
}

proj <- "C:/Users/Axis3/Github Repos/UniWork/IFN646/Project/Task 2"
out_dir <- file.path(proj, "crackling_inputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("Loading gene annotation and genome...\n")
Hs  <- BSgenome.Hsapiens.UCSC.hg38
edb <- EnsDb.Hsapiens.v86

gene_symbol <- "CXCL10"

cat("Finding gene:", gene_symbol, "...\n")
gr_gene <- genes(edb, filter = AnnotationFilter::SymbolFilter(gene_symbol))
stopifnot(length(gr_gene) >= 1)
ensg_id <- gr_gene$gene_id[1]


# Transcripts for the gene
cat("Retrieving transcripts for gene ID:", ensg_id, "...\n")
txs <- transcripts(edb, filter = AnnotationFilter::GeneIdFilter(ensg_id))

## --- CRITICAL HARMONIZATION  ---
cat("Harmonizing transcript coordinates to hg38...\n")
txs <- keepStandardChromosomes(txs, pruning.mode = "coarse")
seqlevelsStyle(txs) <- seqlevelsStyle(Hs)               # make 'chr' names match
seqlevels(txs) <- base::intersect(seqlevels(txs), seqlevels(Hs))
seqinfo(txs)   <- seqinfo(Hs)[seqlevels(txs)]           # attach hg38 lengths

### Sanity check
cat("Transcript summary:\n")
print(summary(txs))


# 1-bp TSS point per transcript (strand aware)
cat("Calculating TSS positions...\n")
tss <- promoters(txs, upstream = 0, downstream = 1)
# Use the actual genomic coordinate of the TSS (start for '+', end for '-')
pos  <- ifelse(as.vector(strand(tss) == "+"), start(tss), end(tss))
plus <- as.vector(strand(tss) == "+")  # <- avoid S4/Rle indexing problem
minus <- !plus

cat("TSS summary:\n")
print(summary(tss))

# Helper to make a window [a,b] around TSS (a,b in bp relative to TSS, strand-aware)
make_tss_window <- function(a, b) {
  n <- length(tss)
  s <- integer(n); e <- integer(n)

  # + strand: TSS + [a,b]
  s[plus]  <- pos[plus] + a
  e[plus]  <- pos[plus] + b
  # - strand: TSS - [b,a]
  s[minus] <- pos[minus] - b
  e[minus] <- pos[minus] - a

  # clamp to chromosome bounds
  chr_len <- as.integer(seqlengths(Hs)[as.character(seqnames(tss))])
  s <- pmax(1L, pmin(s, chr_len))
  e <- pmax(1L, pmin(e, chr_len))

  gr <- GRanges(seqnames = seqnames(tss),
                ranges   = IRanges(start = pmin(s, e), end = pmax(s, e)),
                strand   = strand(tss))
  seqinfo(gr) <- seqinfo(tss)
  trim(gr)  # guarantees width > 0
}

## Examples:
cat("Finding Core and Wide TSS windows...\n")
core <- make_tss_window(25, 150)     # CRISPRi sweet-spot
wide <- make_tss_window(-500, 500)   # cast a wider net for Crackling

# (Optional) deduplicate identical intervals across transcripts
wide <- reduce(wide, ignore.strand = FALSE)

# Write out for Crackling
library(Biostrings)
seqs_core <- getSeq(Hs, core)
names(seqs_core) <- paste0("CXCL10_TSS_core|",
                           as.character(seqnames(core)), ":",
                           start(core), "-", end(core), "(",
                           as.character(strand(core)), ")")
Biostrings::writeXStringSet(seqs_core, file.path(out_dir, "CXCL10_TSS_core.fa"))

seqs_wide <- getSeq(Hs, wide)
names(seqs_wide) <- paste0("CXCL10_TSS_win|",
                           as.character(seqnames(wide)), ":",
                           start(wide), "-", end(wide), "(",
                           as.character(strand(wide)), ")")
Biostrings::writeXStringSet(seqs_wide, file.path(out_dir, "CXCL10_TSS_wide.fa"))


### Sanity checks
core_fp <- file.path(out_dir, "CXCL10_TSS_core.fa")
wide_fp <- file.path(out_dir, "CXCL10_TSS_wide.fa")
print(list(core_exists=file.exists(core_fp), wide_exists=file.exists(wide_fp)))
print(list(core_size=unclass(file.info(core_fp)$size),
           wide_size=unclass(file.info(wide_fp)$size)))


core_fa <- Biostrings::readDNAStringSet(core_fp)
wide_fa <- Biostrings::readDNAStringSet(wide_fp)
cat("core records:", length(core_fa), "wide records:", length(wide_fa), "\n")
cat("core widths  :", paste(unique(Biostrings::width(core_fa)), collapse=", "), "\n")
cat("wide widths  :", paste(unique(Biostrings::width(wide_fa)), collapse=", "), "\n")