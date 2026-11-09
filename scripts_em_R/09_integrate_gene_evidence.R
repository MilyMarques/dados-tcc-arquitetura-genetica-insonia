#!/usr/bin/env Rscript
# Integra MAGMA, mapeamento posicional, eQTL e interação cromatínica por gene e estudo.

options(stringsAsFactors = FALSE, scipen = 999)
arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(arg)) sub("^--file=", "", arg[1]) else "R/09_integrate_gene_evidence.R"
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
in_dir <- file.path(root, "projeto_insonia_tcc", "Resultados", "analise_padronizada")

read_tsv <- function(name) utils::read.delim(file.path(in_dir, name), check.names = FALSE, na.strings = c("", "NA"), quote = "", stringsAsFactors = FALSE)
write_tsv <- function(x, name) utils::write.table(x, file.path(in_dir, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "")
or_false <- function(x) { x[is.na(x)] <- FALSE; as.logical(x) }

magma <- read_tsv("05_magma_all_genes.tsv")
mapped <- read_tsv("04_mapped_genes_and_evidence.tsv")

# Harmoniza nomes e mantém uma linha por gene-estudo para cada fonte de evidência.
m <- magma[, c("GENE", "SYMBOL", "study_id", "label", "evidence_unit", "P", "FDR_BH", "bonferroni_significant", "FDR_significant")]
names(m)[names(m) == "SYMBOL"] <- "symbol_magma"
map <- mapped[, c("ensg", "symbol", "study_id", "label", "evidence_unit", "type", "positional_mapping", "eqtl_mapping", "chromatin_interaction_mapping")]
names(map)[1] <- "GENE"

matrix <- merge(m, map, by = c("GENE", "study_id", "label", "evidence_unit"), all = TRUE, sort = FALSE)
matrix$SYMBOL <- ifelse(!is.na(matrix$symbol_magma) & matrix$symbol_magma != "", matrix$symbol_magma, matrix$symbol)
matrix$symbol_magma <- NULL; matrix$symbol <- NULL
matrix$bonferroni_significant <- or_false(matrix$bonferroni_significant)
matrix$FDR_significant <- or_false(matrix$FDR_significant)
matrix$positional_mapping <- or_false(matrix$positional_mapping)
matrix$eqtl_mapping <- or_false(matrix$eqtl_mapping)
matrix$chromatin_interaction_mapping <- or_false(matrix$chromatin_interaction_mapping)
matrix$mapping_modalities <- matrix$positional_mapping + matrix$eqtl_mapping + matrix$chromatin_interaction_mapping
matrix$any_mapping <- matrix$mapping_modalities > 0
matrix <- matrix[order(matrix$GENE, matrix$study_id), ]

# Resume evidência por gene; Jansen e Watanabe constituem uma unidade UKB única.
by_gene <- split(matrix, matrix$GENE)
summary <- do.call(rbind, lapply(by_gene, function(d) {
  mapped_types <- unique(d$type[!is.na(d$type) & d$type != ""])
  bonf_units <- unique(d$evidence_unit[d$bonferroni_significant])
  fdr_units <- unique(d$evidence_unit[d$FDR_significant])
  pos_units <- unique(d$evidence_unit[d$positional_mapping])
  eqtl_units <- unique(d$evidence_unit[d$eqtl_mapping])
  ci_units <- unique(d$evidence_unit[d$chromatin_interaction_mapping])
  gene_type <- if ("protein_coding" %in% mapped_types) "protein_coding" else if (length(mapped_types)) paste(mapped_types, collapse = ";") else NA_character_
  data.frame(
    GENE = d$GENE[1], SYMBOL = d$SYMBOL[which(!is.na(d$SYMBOL) & d$SYMBOL != "")[1]],
    gene_type = gene_type, protein_coding = identical(gene_type, "protein_coding"),
    n_fuma_studies_magma_bonf = length(unique(d$study_id[d$bonferroni_significant])),
    n_independent_units_magma_bonf = length(bonf_units), magma_bonf_units = paste(sort(bonf_units), collapse = ";"),
    n_independent_units_magma_fdr = length(fdr_units),
    positional = length(pos_units) > 0, positional_units = paste(sort(pos_units), collapse = ";"),
    eqtl_brain = length(eqtl_units) > 0, eqtl_units = paste(sort(eqtl_units), collapse = ";"),
    chromatin_cortex = length(ci_units) > 0, chromatin_units = paste(sort(ci_units), collapse = ";"),
    n_mapping_modalities = sum(length(pos_units) > 0, length(eqtl_units) > 0, length(ci_units) > 0),
    n_independent_units_with_any_mapping = length(unique(d$evidence_unit[d$any_mapping])),
    best_magma_p = if (any(is.finite(d$P))) min(d$P, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}))

# Classificação pré-especificada: o MAGMA replicado é o núcleo; mapeamentos dão suporte funcional.
summary$priority_tier <- "exploratory"
replicated <- summary$protein_coding & summary$n_independent_units_magma_bonf >= 2
summary$priority_tier[replicated & summary$n_mapping_modalities >= 2] <- "high_priority"
summary$priority_tier[replicated & summary$n_mapping_modalities == 1] <- "intermediate_priority"
summary$priority_tier[replicated & summary$n_mapping_modalities == 0] <- "replicated_without_mapping"
summary$priority_tier[summary$protein_coding & summary$n_independent_units_magma_bonf == 1 & summary$n_mapping_modalities >= 2] <- "exploratory_functional"
summary <- summary[order(factor(summary$priority_tier, levels = c("high_priority", "intermediate_priority", "replicated_without_mapping", "exploratory_functional", "exploratory")), -summary$n_independent_units_magma_bonf, summary$best_magma_p), ]

follow_up <- summary[summary$priority_tier %in% c("high_priority", "intermediate_priority", "replicated_without_mapping"), ]
audit <- data.frame(
  metric = c("gene_study_rows", "genes_total", "protein_coding_genes", "high_priority", "intermediate_priority", "replicated_without_mapping"),
  value = c(nrow(matrix), nrow(summary), sum(summary$protein_coding), sum(summary$priority_tier == "high_priority"), sum(summary$priority_tier == "intermediate_priority"), sum(summary$priority_tier == "replicated_without_mapping"))
)

write_tsv(matrix, "09_gene_x_evidence_x_study.tsv")
write_tsv(summary, "10_gene_evidence_summary.tsv")
write_tsv(follow_up, "11_provisional_consensus_protein_coding.tsv")
write_tsv(audit, "12_consensus_audit.tsv")
cat("Matriz e consenso provisório gerados em: ", in_dir, "\n", sep = "")
print(audit)
