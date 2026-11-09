#!/usr/bin/env Rscript
# Documenta a transição entre a análise anterior Jansen/Watanabe e a análise FUMA atual.

options(stringsAsFactors = FALSE)
arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(arg)) sub("^--file=", "", arg[1]) else "R/11_compare_legacy_jansen_watanabe.R"
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
legacy_dir <- file.path(root, "projeto_insonia_tcc")
out_dir <- file.path(legacy_dir, "Resultados", "analise_padronizada")
old_j <- unique(utils::read.delim(file.path(legacy_dir, "magma_genes_significativos_jansen.tsv"), stringsAsFactors = FALSE)$SYMBOL)
old_w <- unique(utils::read.delim(file.path(legacy_dir, "magma_genes_significativos_watanabe.tsv"), stringsAsFactors = FALSE)$SYMBOL)
new <- utils::read.delim(file.path(out_dir, "06_magma_bonferroni_genes.tsv"), stringsAsFactors = FALSE)
new_j <- unique(new$SYMBOL[new$study_id == "jansen_2019_public"])
new_w <- unique(new$SYMBOL[new$study_id == "watanabe_2022_ukb"])

old_union <- union(old_j, old_w); new_union <- union(new_j, new_w)
metrics <- data.frame(
  metric = c("legacy_jansen", "legacy_watanabe", "legacy_intersection", "legacy_union", "current_jansen", "current_watanabe", "current_intersection", "current_union", "legacy_current_union_intersection"),
  n_genes = c(length(old_j), length(old_w), length(intersect(old_j, old_w)), length(old_union), length(new_j), length(new_w), length(intersect(new_j, new_w)), length(new_union), length(intersect(old_union, new_union)))
)
all_genes <- sort(union(old_union, new_union))
membership <- data.frame(
  SYMBOL = all_genes,
  legacy_jansen = all_genes %in% old_j, legacy_watanabe = all_genes %in% old_w,
  current_jansen = all_genes %in% new_j, current_watanabe = all_genes %in% new_w,
  legacy_union = all_genes %in% old_union, current_union = all_genes %in% new_union,
  transition = ifelse(all_genes %in% old_union & all_genes %in% new_union, "retained", ifelse(all_genes %in% old_union, "legacy_only", "current_only")),
  stringsAsFactors = FALSE
)
utils::write.table(metrics, file.path(out_dir, "17_legacy_vs_current_jw_metrics.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
utils::write.table(membership, file.path(out_dir, "18_legacy_vs_current_jw_gene_membership.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(metrics)
