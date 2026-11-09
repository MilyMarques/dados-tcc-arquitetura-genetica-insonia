#!/usr/bin/env Rscript
# Integra o núcleo de 16 genes com evidência SNP2GENE, MAGMA, sensibilidade
# OpenGWAS-UKB e resultados FDR-significativos do GENE2FUNC.

options(stringsAsFactors = FALSE, scipen = 999)
arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(arg)) sub("^--file=", "", arg[1]) else "R/14_prepare_final_gene_evidence_table.R"
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(root, "projeto_insonia_tcc", "Resultados", "analise_padronizada")
g2f_dir <- file.path(root, "projeto_insonia_tcc", "Resultados_GENE2FUNC", "FUMA_gene2func761816_16_genes")
enrichment_dir <- file.path(root, "projeto_insonia_tcc", "enrichment_inputs")
dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                    na.strings = c("", "NA"), quote = "")
}
write_tsv <- function(x, name) {
  utils::write.table(x, file.path(analysis_dir, name), sep = "\t", row.names = FALSE,
                     quote = FALSE, na = "")
}
collapse_or_blank <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (!length(x)) "" else paste(x, collapse = "; ")
}

core <- read_tsv(file.path(analysis_dir, "11_provisional_consensus_protein_coding.tsv"))
ukb <- read_tsv(file.path(analysis_dir, "21_core16_opengwas_ukb_sensitivity_summary.tsv"))
gs <- read_tsv(file.path(g2f_dir, "GS.txt"))
gs$adjP <- suppressWarnings(as.numeric(gs$adjP))
sig_gs <- gs[is.finite(gs$adjP) & gs$adjP < 0.05, ]

# Associa cada gene aos conjuntos significativos. A participação não equivale a
# causalidade; ela apenas documenta quais genes dirigem cada resultado de enrichment.
g2f_by_symbol <- lapply(core$SYMBOL, function(symbol) {
  # Campos de genes vazios/NA não devem gerar uma linha artificial na seleção.
  member <- vapply(sig_gs$genes, function(x) {
    if (is.na(x) || x == "") return(FALSE)
    symbol %in% strsplit(x, ":", fixed = TRUE)[[1]]
  }, logical(1))
  hits <- sig_gs[member, , drop = FALSE]
  data.frame(
    SYMBOL = symbol,
    gene2func_n_fdr_gene_sets = nrow(hits),
    gene2func_fdr_gene_sets = if (nrow(hits)) collapse_or_blank(paste0(hits$Category, ":", hits$GeneSet)) else "",
    gene2func_min_fdr = if (nrow(hits)) min(hits$adjP) else NA_real_,
    stringsAsFactors = FALSE
  )
})
g2f_by_symbol <- do.call(rbind, g2f_by_symbol)

final <- merge(core, ukb[, c("GENE", "n_ukb_sensitivity_runs_tested",
                              "n_ukb_sensitivity_runs_bonferroni",
                              "n_ukb_sensitivity_runs_fdr", "ukb_runs_bonferroni",
                              "best_ukb_magma_p")], by = "GENE", all.x = TRUE)
final <- merge(final, g2f_by_symbol, by = "SYMBOL", all.x = TRUE)
final$main_evidence_summary <- paste0(
  "MAGMA Bonferroni em ", final$n_independent_units_magma_bonf,
  " unidade(s) independente(s): ", final$magma_bonf_units,
  "; modalidades SNP2GENE: ", final$n_mapping_modalities
)
final$ukb_sensitivity_summary <- paste0(
  "FDR em ", final$n_ukb_sensitivity_runs_fdr, "/", final$n_ukb_sensitivity_runs_tested,
  "; Bonferroni em ", final$n_ukb_sensitivity_runs_bonferroni, "/", final$n_ukb_sensitivity_runs_tested,
  " (sensibilidade UKB, não independente)"
)

priority_order <- c("high_priority", "intermediate_priority")
final <- final[order(match(final$priority_tier, priority_order), -final$n_independent_units_magma_bonf,
                     final$best_magma_p), ]

keep <- c("GENE", "SYMBOL", "priority_tier", "main_evidence_summary",
          "n_independent_units_magma_bonf", "magma_bonf_units", "best_magma_p",
          "positional", "eqtl_brain", "chromatin_cortex", "n_mapping_modalities",
          "ukb_sensitivity_summary", "best_ukb_magma_p",
          "gene2func_n_fdr_gene_sets", "gene2func_min_fdr", "gene2func_fdr_gene_sets")
write_tsv(final[, keep], "22_core16_final_evidence_table.tsv")

g2f_summary <- sig_gs[, c("Category", "GeneSet", "N_genes", "N_overlap", "p", "adjP", "genes")]
g2f_summary <- g2f_summary[order(g2f_summary$adjP), ]
write_tsv(g2f_summary, "23_gene2func_fdr_significant_gene_sets.tsv")

meta <- data.frame(
  item = c("core_definition", "main_studies", "independence_rule", "SNP2GENE_evidence",
           "UKB_role", "GENE2FUNC_input", "GENE2FUNC_background", "GENE2FUNC_interpretation"),
  value = c(
    "16 protein-coding genes with MAGMA Bonferroni significance in >=2 independent units and >=1 SNP2GENE mapping modality",
    "Hammerschlag 2017; Lane 2019; Jansen 2019 public UKB; Watanabe 2022 public UKB",
    "Jansen and Watanabe are counted as one UKB unit because of substantial sample overlap",
    "positive positional, brain eQTL and/or cortical chromatin-interaction mapping only",
    "ukb-a-13 and ukb-b-3957 are sensitivity analyses and do not change the core",
    "16 Ensembl gene IDs",
    "18,717 genes testable by MAGMA in all four main studies",
    "membership in a significant gene set identifies contributors to enrichment, not causal genes"
  ), stringsAsFactors = FALSE
)
write_tsv(meta, "24_final_evidence_table_readme.tsv")

# Listas simples, sem cabeçalho: adequadas para colar ou carregar em g:Profiler e ShinyGO.
# Símbolos HGNC são os mais convenientes nas interfaces; IDs Ensembl são fornecidos
# para uma alternativa inequívoca de identificação.
background <- read_tsv(file.path(root, "projeto_insonia_tcc", "gene2func_inputs",
                                 "gene2func_magma_testable_background_ensembl.tsv"))
magma <- read_tsv(file.path(analysis_dir, "05_magma_all_genes.tsv"))
symbol_lookup <- magma[match(background$GENE, magma$GENE), c("GENE", "SYMBOL")]
symbol_lookup <- symbol_lookup[!duplicated(symbol_lookup$GENE), ]
writeLines(core$SYMBOL, file.path(enrichment_dir, "core16_gene_symbols.txt"), useBytes = TRUE)
writeLines(core$GENE, file.path(enrichment_dir, "core16_ensembl_ids.txt"), useBytes = TRUE)
writeLines(background$GENE, file.path(enrichment_dir, "magma_background_18717_ensembl_ids.txt"), useBytes = TRUE)
writeLines(symbol_lookup$SYMBOL[!is.na(symbol_lookup$SYMBOL) & symbol_lookup$SYMBOL != ""],
           file.path(enrichment_dir, "magma_background_18717_gene_symbols.txt"), useBytes = TRUE)

cat("Tabelas finais geradas em: ", analysis_dir, "\n", sep = "")
cat("- 22_core16_final_evidence_table.tsv\n- 23_gene2func_fdr_significant_gene_sets.tsv\n- 24_final_evidence_table_readme.tsv\n")
cat("Listas para enrichment geradas em: ", enrichment_dir, "\n", sep = "")
