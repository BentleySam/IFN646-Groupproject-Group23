# ---- crackling_summarise.R -----------------------------
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
})

# >>> EDIT THIS <<<  (or pass via commandArgs if you prefer)
in_file <- "C:/Users/Axis3/Github Repos/UniWork/IFN646/Project/Task 2/cracklingOutput/Task_2_CXCL10_relaxed_avg_65_wide-guides.txt"

# --------------------------------------------------------
# Read & clean
df <- readr::read_csv(in_file, na = c("?", "NA", ""), show_col_types = FALSE)

# Parse the header like: CXCL10_TSS_win|chr4:76022997-76023997(-)
hdr <- df |>
  tidyr::extract(
    header,
    into = c("region_id","hdr_chr","hdr_start","hdr_end","hdr_strand"),
    regex = "^([^|]+)\\|(chr[^:]+):(\\d+)-(\\d+)\\(([+-])\\)$",
    remove = FALSE
  ) |>
  mutate(
    hdr_start = as.integer(hdr_start),
    hdr_end   = as.integer(hdr_end)
  )

# Derive spacer & PAM (Crackling’s `seq` ends with PAM; usually NGG for SpCas9)
hdr <- hdr |>
  mutate(
    seq_len   = nchar(seq),
    pam       = substr(seq, seq_len - 2, seq_len),
    spacer20  = substr(seq, 1, pmax(0, seq_len - 3))
  )

# Absolute genomic coords for the protospacer (based on header region + local start/end)
# If the REGION is '-' we flip coordinates relative to the region end.
hdr <- hdr |>
  mutate(
    abs_start = ifelse(hdr_strand == "+", hdr_start + start - 1L, hdr_end - end   + 1L),
    abs_end   = ifelse(hdr_strand == "+", hdr_start + end   - 1L, hdr_end - start + 1L),
    protospacer_strand = ifelse(hdr_strand == "+", strand, ifelse(strand == "+", "-", "+")),
    region_center = floor((hdr_start + hdr_end)/2),
    dist_to_center = pmin(abs(abs_start - region_center), abs(abs_end - region_center))
  )

# Turn 0/1 into logicals for the main pass/fail flags
bin_cols <- c("isUnique","passedG20","passedTTTT","passedATPercent","passedSecondaryStructure",
              "acceptedByMm10db","acceptedBySgRnaScorer","passedBowtie","passedOffTargetScore",
              "passedAvoidLeadingT")
hdr <- hdr |>
  mutate(across(any_of(bin_cols), ~ ifelse(is.na(.x), NA, .x == 1)))

# Build a compact reason string for anything that fails
reason_cols <- c(
  G20 = "passedG20", TTTT = "passedTTTT", ATpct = "passedATPercent",
  RNAfold = "passedSecondaryStructure", mm10db = "acceptedByMm10db",
  sgRNAScorer = "acceptedBySgRnaScorer", Bowtie = "passedBowtie",
  OffTarget = "passedOffTargetScore", leadT = "passedAvoidLeadingT"
)

fail_reasons <- function(row) {
  out <- character()
  for (nm in names(reason_cols)) {
    col <- reason_cols[[nm]]
    val <- row[[col]]
    if (isFALSE(val)) out <- c(out, nm)
  }
  if (length(out) == 0) NA_character_ else paste(out, collapse = "|")
}

hdr$fail_reason <- pmap_chr(as.list(hdr), fail_reasons)

# A single “keep” flag matching your AVG-65 run idea
# (accepted by both efficiency models; consensus >=2; Bowtie & ISSL passed)
hdr <- hdr |>
  mutate(
    consensus_ok = !is.na(consensusCount) & consensusCount >= 2,
    efficiency_ok = acceptedByMm10db & acceptedBySgRnaScorer,
    bowtie_ok = passedBowtie,
    offtarget_ok = passedOffTargetScore,
    keep = consensus_ok & efficiency_ok & bowtie_ok & offtarget_ok
  )

# Ranking: prefer keep==TRUE, then higher cfd/mit/sgRNAScorer2
hdr <- hdr |>
  mutate(
    cfd_num = as.numeric(cfdOfftargetscore),
    mit_num = as.numeric(mitOfftargetscore),
    sgrna2_num = suppressWarnings(as.numeric(sgrnascorer2score)),
    rank_key = dplyr::dense_rank(desc(keep)) * 1e9 +
               dplyr::dense_rank(desc(replace_na(cfd_num, -Inf))) * 1e6 +
               dplyr::dense_rank(desc(replace_na(mit_num, -Inf))) * 1e3 +
               dplyr::dense_rank(desc(replace_na(sgrna2_num, -Inf)))
  ) |>
  arrange(desc(keep), desc(cfd_num), desc(mit_num), desc(sgrna2_num), dist_to_center)

# Select tidy columns
tidy <- hdr |>
  transmute(
    spacer = spacer20, pam, seq,
    genome_chr = bowtieChr %||% hdr_chr,   # prefer Bowtie chr if present
    abs_start, abs_end, protospacer_strand, region = header,
    dist_to_center,
    isUnique, consensusCount, efficiency_ok, bowtie_ok, offtarget_ok, keep,
    sgrnascorer2 = sgrna2_num,
    mit = mit_num,
    cfd = cfd_num,
    fail_reason
  )

# Write outputs next to the input
out_dir <- dirname(in_file)
stem <- sub("\\.txt$", "", basename(in_file))
summary_csv <- file.path(out_dir, paste0(stem, "_summary_clean.csv"))
top_csv     <- file.path(out_dir, paste0(stem, "_TOP10.csv"))
bed_file    <- file.path(out_dir, paste0(stem, ".bed"))

readr::write_csv(tidy, summary_csv)

# Top 10 keep==TRUE (else top 10 overall)
top <- tidy |> filter(keep) |> head(10)
if (nrow(top) == 0) top <- tidy |> head(10)
readr::write_csv(top, top_csv)

# Optional BED (0-based, half-open). Uses abs_start/abs_end as-is and bowtieChr if available.
# Score: scaled CFD in [0..1000].
bed <- hdr |>
  mutate(
    bed_chr  = bowtieChr %||% hdr_chr,
    bed_start0 = pmax(0L, abs_start - 1L),
    bed_end0   = abs_end,
    bed_name = paste0(spacer20, "_", pam, "_", protospacer_strand),
    bed_score = pmax(0, pmin(1000, round(replace_na(cfd_num, 0) / 100 * 1000))),
    bed_strand = ifelse(protospacer_strand == "+", "+", "-")
  ) |>
  select(bed_chr, bed_start0, bed_end0, bed_name, bed_score, bed_strand)

readr::write_tsv(bed, bed_file, col_names = FALSE)

message("Wrote:\n  ", summary_csv, "\n  ", top_csv, "\n  ", bed_file)
# --------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || all(is.na(a))) b else a
# ---- end ------------------------------------------------
