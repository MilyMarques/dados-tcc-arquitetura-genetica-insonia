#!/usr/bin/env Rscript
# Prepara lista consenso e universo testável MAGMA para GENE2FUNC e outras ferramentas.

options(stringsAsFactors = FALSE)
arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(arg)) sub("^--file=", "", arg[1]) else "R/13_prepare_gene2func_inputs.R"
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
source_dir <- file.path(root, "projeto_insonia_tcc", "Resultados", "analise_padronizada")
out_dir <- file.path(root, "projeto_insonia_tcc", "gene2func_inputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

core <- utils::read.delim(file.path(source_dir, "11_provisional_consensus_protein_coding.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
magma <- utils::read.delim(file.path(source_dir, "05_magma_all_genes.tsv"), check.names = FALSE, stringsAsFactors = FALSE)

# Universo: genes MAGMA efetivamente testados em todos os quatro jobs principais.
# Isso garante que cada gene do background teria tido oportunidade de entrar no consenso.
primary_studies <- c("hammerschlag_2017", "lane_2019_extreme", "jansen_2019_public", "watanabe_2022_ukb")
tested_sets <- lapply(primary_studies, function(study) unique(magma$GENE[magma$study_id == study & is.finite(magma$P)]))
background_ids <- Reduce(intersect, tested_sets)
symbol_map <- magma[match(background_ids, magma$GENE), c("GENE", "SYMBOL")]
symbol_map <- symbol_map[!duplicated(symbol_map$GENE), ]
background_symbols <- symbol_map$SYMBOL[!is.na(symbol_map$SYMBOL) & symbol_map$SYMBOL != ""]

write_one_column <- function(values, file) {
  utils::write.table(data.frame(GENE = values), file, sep = "\t", row.names = FALSE, quote = FALSE)
}
write_one_column(core$GENE, file.path(out_dir, "gene2func_core16_ensembl.tsv"))
write_one_column(core$SYMBOL, file.path(out_dir, "gene2func_core16_symbols.tsv"))
write_one_column(background_ids, file.path(out_dir, "gene2func_magma_testable_background_ensembl.tsv"))
write_one_column(background_symbols, file.path(out_dir, "gene2func_magma_testable_background_symbols.tsv"))

manifest <- data.frame(
  item = c("analysis_role", "input_list", "input_gene_count", "input_selection_rule", "background_list", "background_gene_count", "background_definition", "identifier_recommendation", "important_note"),
  value = c(
    "GENE2FUNC primary exploratory functional interpretation",
    "gene2func_core16_ensembl.tsv",
    nrow(core),
    "protein-coding genes, MAGMA Bonferroni significant in >=2 independent units; Jansen and Watanabe count as one UKB unit; >=1 mapping modality",
    "gene2func_magma_testable_background_ensembl.tsv",
    length(background_ids),
    "intersection of genes with finite MAGMA P-values in Hammerschlag, Lane, Jansen and Watanabe; represents the gene universe eligible for selection",
    "Ensembl IDs are primary; symbol files supplied only for interfaces that require symbols",
    "The core has 16 protein-coding genes. The custom background is MAGMA-testable universe, which is statistically preferable to an all-genome background."
  ), stringsAsFactors = FALSE
)
utils::write.table(manifest, file.path(out_dir, "gene2func_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("Arquivos GENE2FUNC preparados em: ", out_dir, "\n", sep = "")
print(manifest)
