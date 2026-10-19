suppressPackageStartupMessages(library(DIALOGUE))
D <- "/LARGE1/gr10634/gaozy/aml_niche_net/ws_dialogue"
pick <- function(dir) { f <- list.files(dir, pattern = "^DLG\\.full\\.output.*\\.rds$", full.names = TRUE); if (!length(f)) f <- list.files(dir, pattern = "^DIALOGUE2_.*\\.rds$", full.names = TRUE); f[1] }
smp <- function(R) { if (!is.null(R$samples)) sort(unique(as.character(R$samples))) else sort(Reduce(union, lapply(R$scores, function(S) unique(as.character(S$samples))))) }
cat("=== stability draws: sample sets and emp.p fingerprints ===\n")
for (run in c("R1", "R2")) { prev <- NULL
  for (b in 1:20) { f <- pick(file.path(D, "stability", paste0(run, "_stab", b))); if (is.na(f)) { cat(run, b, "no file\n"); next }
    R <- readRDS(f); s <- smp(R); ep <- R$emp.p; fp <- if (is.null(ep)) "NULL" else paste(round(as.numeric(as.matrix(ep)), 3), collapse = ",")
    cat(sprintf("%s stab%02d  n_samples %2d  same-set-as-previous %-5s  emp.p [%s]  file %s\n", run, b, length(s), if (is.null(prev)) "-" else identical(s, prev), substr(fp, 1, 60), basename(f))); prev <- s } }
cat("\n=== shuffles: emp.p fingerprints and identity with shuffle 1 ===\n"); E1 <- NULL
for (b in 1:20) { f <- pick(file.path(D, "shuffle", paste0("R1_shuf", b))); if (is.na(f)) { cat("shuf", b, "no file\n"); next }
  R <- readRDS(f); ep <- as.matrix(R$emp.p); if (b == 1) E1 <- ep
  cat(sprintf("shuf%02d  emp.p [%s]  identical-to-shuf1 %s  MCP.cell.types %s\n", b, paste(round(as.numeric(ep), 3), collapse = ","), identical(ep, E1), paste(sapply(R$MCP.cell.types, function(x) paste(x, collapse = "+")), collapse = " | "))) }
