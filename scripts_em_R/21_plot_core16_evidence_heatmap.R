# Heatmap reproducível das evidências que definem o Core16.
# Entrada: projeto_insonia_tcc/Resultados/analise_padronizada/22_core16_final_evidence_table.tsv
# Saída: projeto_insonia_tcc/Figuras_TCC/core16_evidence_heatmap.png
# Usa apenas R base.

root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
input_file <- file.path(root, "projeto_insonia_tcc", "Resultados", "analise_padronizada", "22_core16_final_evidence_table.tsv")
output_dir <- file.path(root, "projeto_insonia_tcc", "Figuras_TCC")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "core16_evidence_heatmap.png")

d <- read.delim(input_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("SYMBOL", "n_independent_units_magma_bonf", "positional", "eqtl_brain", "chromatin_cortex") %in% names(d)))

# A primeira coluna representa 2 ou 3 unidades independentes com MAGMA-Bonferroni.
# As demais representam presença/ausência de evidência SNP2GENE usada na regra do Core16.
mat <- cbind(
  `MAGMA\n(unidades)` = d$n_independent_units_magma_bonf,
  Posicional = as.integer(d$positional == TRUE),
  `eQTL\ncerebral` = as.integer(d$eqtl_brain == TRUE),
  `Cromatina\ncortical` = as.integer(d$chromatin_cortex == TRUE)
)
rownames(mat) <- d$SYMBOL

png(output_file, width = 2400, height = 1900, res = 300)
par(mar = c(3.5, 8.5, 4.5, 1), family = "sans")
n <- nrow(mat); k <- ncol(mat)
plot(NA, xlim = c(0, k), ylim = c(0, n), axes = FALSE, xlab = "", ylab = "", bty = "n")
for (i in seq_len(n)) {
  for (j in seq_len(k)) {
    value <- mat[i, j]
    col <- if (j == 1 && value == 3) "#1F4E79" else if (j == 1 && value == 2) "#5B9BD5" else if (j > 1 && value == 1) "#5B9BD5" else "#F1F1F1"
    rect(j - 1, n - i, j, n - i + 1, col = col, border = "white", lwd = 2)
    label <- if (j == 1) as.character(value) else if (value == 1) "Sim" else "Não"
    text(j - .5, n - i + .5, label, cex = .85, col = if (col == "#1F4E79") "white" else "black")
  }
}
axis(2, at = n - seq_len(n) + .5, labels = rownames(mat), las = 1, tick = FALSE, cex.axis = .9)
axis(3, at = seq_len(k) - .5, labels = colnames(mat), tick = FALSE, cex.axis = .92)
title("Evidências integradas que definem o Core16", cex.main = 1.35, line = 2.1)
mtext("MAGMA: número de unidades independentes com significância por Bonferroni; demais colunas: modalidades SNP2GENE.", side = 1, line = 1.6, cex = .72)
dev.off()

message("Heatmap salvo em: ", output_file)
