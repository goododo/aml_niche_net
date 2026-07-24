#!/usr/bin/env Rscript
# 01_infercnv_make_gene_order.R ----
# Build the inferCNV gene_order_file from the 10x GRCh38-2020-A genes.gtf.
#
# Pure-R decompress + filter (no gunzip/awk/zcat/shell): the 10x GTF is BGZF (block-gzip), which
# fread's gzip reader and some `zcat` builds can't handle. R's gzfile() connection (zlib) reads
# BGZF's concatenated blocks fine, so we read in chunks and keep only `gene` feature lines.
#
# INPUT  : args[1] = /path/to/refdata-gex-GRCh38-2020-A/genes/genes.gtf
# OUTPUT : args[2] = gene_order file (no header, ordered by chromosome then start):
#            gene <TAB> chr <TAB> start <TAB> end
#          This is the INFERCNV_GENE_ORDER consumed by 40/44.
#
# 3a-1 rewire (no logic change): here-anchored config source; bare fwrite -> fwrite_safe. The
#   positional-arg CLI is unchanged. (config-completion candidate: once INFERCNV_GENE_ORDER lands
#   in config_malignancy.R, args[2] could default to it; kept required for now.)
#
# Run (from PROJECT ROOT so here() resolves via .here):
#   Rscript scripts/02_malignancy/01_infercnv_make_gene_order.R \
#     /path/to/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
#     /LARGE1/gr10634/gaozy/aml_niche_net/reference/gencode_GRCh38_gene_order.txt

suppressPackageStartupMessages(library(data.table))

## ---- config (here-anchored; zero hard-coded paths) ----
suppressPackageStartupMessages({ library(here) })
source(here::here("scripts", "config", "config_paths.R"))
source(here::here("scripts", "config", "config_malignancy.R"))
source(here::here("scripts", "config", "utils.R"))

args <- commandArgs(TRUE)
gtf  <- args[1]; out <- args[2]
stopifnot(length(args) == 2, file.exists(gtf))

# detect gzip by magic bytes (BGZF or plain gzip both start 1f 8b)
con0  <- file(gtf, "rb"); magic <- readBin(con0, "raw", 2L); close(con0)
is_gz <- length(magic) == 2L && magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b)
message(sprintf("[gtf] %s -> reading as %s", basename(gtf), if (is_gz) "gzip/BGZF" else "plain text"))

# stream the (decompressed) GTF in chunks; keep only `gene` feature lines.
con <- if (is_gz) gzfile(gtf, "rt") else file(gtf, "rt")
CHUNK <- 2e6L
gene_lines <- character(0)
n_read <- 0L
repeat {
  ll <- readLines(con, n = CHUNK, warn = FALSE)
  if (!length(ll)) break
  n_read <- n_read + length(ll)
  hit <- grepl("\tgene\t", ll, fixed = TRUE)   # feature column == gene (tab-gene-tab)
  if (any(hit)) gene_lines <- c(gene_lines, ll[hit])
}
close(con)
message(sprintf("[gtf] scanned %d lines, found %d gene rows", n_read, length(gene_lines)))
if (!length(gene_lines))
  stop("no 'gene' feature lines found -- is this really a GTF? (check the gzfile readLines test)")

# parse the small set of gene lines
g <- fread(text = gene_lines, sep = "\t", header = FALSE, quote = "",
           col.names = c("chr","src","feature","start","end","score","strand","frame","attr"))
g[, gene := sub('.*gene_name "([^"]+)".*', "\\1", attr)]

# keep main chromosomes (handles "1" and "chr1" naming)
main <- c(as.character(1:22), "X", "Y")
g[, chr_bare := sub("^chr", "", as.character(chr))]
g <- g[chr_bare %in% main]

# one row per gene symbol (first occurrence); order by chromosome + start
g <- g[!duplicated(gene)]
g[, chr_ord := match(chr_bare, main)]
setorder(g, chr_ord, start)

fwrite_safe(g[, .(gene, chr, start, end)], out, sep = "\t", col.names = FALSE)
cat(sprintf("wrote %d genes -> %s\n", nrow(g), out))
cat("first rows:\n"); print(head(g[, .(gene, chr, start, end)]))
