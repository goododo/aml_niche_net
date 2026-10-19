suppressPackageStartupMessages(library(DIALOGUE))
ns <- asNamespace("DIALOGUE"); cat("######## generic.vector2mat ########\n"); print(get("generic.vector2mat", ns))
D <- "/LARGE1/gr10634/gaozy/aml_niche_net/ws_dialogue"
ict_fixed <- function(R) {   # same algorithm, but the logical index is applied in p1's own row order
  emp.p <- R$emp.p; emp.p2 <- R$emp.p
  colnames(emp.p2) <- paste(DIALOGUE:::get.strsplit(colnames(emp.p), "_", 2), DIALOGUE:::get.strsplit(colnames(emp.p), "_", 1), sep = "_")
  emp.p <- cbind(emp.p, emp.p2); out <- list()
  for (x in rownames(emp.p)) { p1 <- DIALOGUE:::generic.vector2mat(emp.p[x, ]); p1 <- p1[R$cell.types, R$cell.types]; diag(p1) <- 0.05
    cell.types <- R$cell.types; lb <- min(rowSums(p1 < 0.1))
    while (lb < ceiling(length(cell.types) / 2)) { drop <- rownames(p1)[rowSums(p1 < 0.1) <= lb][1]; cell.types <- setdiff(cell.types, drop); p1 <- p1[cell.types, cell.types, drop = FALSE]; lb <- min(rowSums(p1 < 0.1)) }
    out[[x]] <- cell.types }
  out
}
for (run in c("R1", "R2", "planted")) { R <- readRDS(file.path(D, paste0("R_", run, ".rds")))
  p1 <- DIALOGUE:::generic.vector2mat(cbind(R$emp.p, setNames(R$emp.p, NULL))[1, ])
  cat("\n=====", run, "| cell.types order:", paste(R$cell.types, collapse = ","), "| generic.vector2mat row order:", paste(rownames(DIALOGUE:::generic.vector2mat(R$emp.p[1, ])), collapse = ","), "\n")
  a <- DIALOGUE:::DIALOGUE.identify.cell.types(R); b <- ict_fixed(R)
  for (m in names(a)) cat(sprintf("  %s  DIALOGUE: %-14s  order-fixed: %-14s  %s\n", m, paste(a[[m]], collapse = "+"), paste(b[[m]], collapse = "+"), if (identical(sort(a[[m]]), sort(b[[m]]))) "" else "<-- DIFFERS")) }
