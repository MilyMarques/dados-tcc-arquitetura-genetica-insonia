#!/usr/bin/env Rscript
# Prepara a análise de sensibilidade histórica baseada nas rodagens preliminares
# de Jansen e Watanabe: união de genes MAGMA significativos e background compatível.

options(stringsAsFactors = FALSE, scipen = 999)
arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(arg)) sub("^--file=", "", arg[1]) else "R/15_prepare_legacy_jansen_watanabe_sensitivity_inputs.R"
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
base <- file.path(root, "projeto_insonia_tcc")
out <- file.path(base, "enrichment_inputs", "legacy_jansen_watanabe_28")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                                              na.strings = c("", "NA"), quote = "")
write_one <- function(x, path) writeLines(as.character(x), path, useBytes = TRUE)

j_sig <- read_tsv(file.path(base, "magma_genes_significativos_jansen.tsv"))
w_sig <- read_tsv(file.path(base, "magma_genes_significativos_watanabe.tsv"))
j_all <- read_tsv(file.path(base, "fuma_results_jansen", "SNP2GENE - Jansen", "jansen_magma.genes.out"))
w_all <- read_tsv(file.path(base, "fuma_results_watanabe", "SNP2GENE - Watanabe", "watanabe_magma.genes.out"))

# A lista histórica é a união dos 26 genes Bonferroni-significativos em cada
# rodada preliminar: 24 compartilhados + 2 exclusivos em cada estudo = 28.
legacy <- rbind(j_sig[, c("GENE", "SYMBOL", "P")], w_sig[, c("GENE", "SYMBOL", "P")])
legacy$P <- suppressWarnings(as.numeric(legacy$P))
legacy <- legacy[order(legacy$GENE, legacy$P), ]
legacy <- legacy[!duplicated(legacy$GENE), ]

# Universo de seleção histórico: genes com p MAGMA finito nas duas rodagens
# preliminares. Usar a interseção evita comparar os 28 genes contra um universo
# diferente do que efetivamente poderia tê-los produzido.
j_all$P <- suppressWarnings(as.numeric(j_all$P))
w_all$P <- suppressWarnings(as.numeric(w_all$P))
j_test <- j_all[is.finite(j_all$P), c("GENE", "SYMBOL")]
w_test <- w_all[is.finite(w_all$P), c("GENE", "SYMBOL")]
bg_ids <- intersect(j_test$GENE, w_test$GENE)
bg_symbols <- j_test$SYMBOL[match(bg_ids, j_test$GENE)]

write_one(legacy$GENE, file.path(out, "legacy_jansen_watanabe_28_ensembl_ids.txt"))
write_one(legacy$SYMBOL, file.path(out, "legacy_jansen_watanabe_28_gene_symbols.txt"))
write_one(bg_ids, file.path(out, "legacy_jansen_watanabe_magma_background_ensembl_ids.txt"))
write_one(bg_symbols[!is.na(bg_symbols) & bg_symbols != ""],
          file.path(out, "legacy_jansen_watanabe_magma_background_gene_symbols.txt"))

membership <- merge(
  merge(legacy[, c("GENE", "SYMBOL", "P")], j_sig[, c("GENE", "P")], by = "GENE", all.x = TRUE, suffixes = c("_best", "_jansen")),
  w_sig[, c("GENE", "P")], by = "GENE", all.x = TRUE
)
membership$in_jansen <- is.finite(suppressWarnings(as.numeric(membership$P_jansen)))
membership$in_watanabe <- is.finite(suppressWarnings(as.numeric(membership$P)))
membership$membership <- ifelse(membership$in_jansen & membership$in_watanabe, "shared",
                                ifelse(membership$in_jansen, "Jansen_only", "Watanabe_only"))
membership$P_watanabe <- membership$P
membership$P <- NULL
membership <- membership[order(membership$membership, membership$SYMBOL), ]
utils::write.table(membership, file.path(out, "legacy_jansen_watanabe_28_membership.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE, na = "")

manifest <- data.frame(
  item = c("analysis_role", "input_definition", "input_gene_count", "shared_genes", "Jansen_only", "Watanabe_only",
           "background_definition", "background_gene_count", "independence_limitation", "recommended_platform_settings"),
  value = c(
    "Historical sensitivity analysis; does not redefine the current 16-gene core",
    "Union of Bonferroni-significant MAGMA genes from preliminary Jansen and Watanabe FUMA jobs",
    nrow(legacy),
    sum(membership$membership == "shared"), sum(membership$membership == "Jansen_only"), sum(membership$membership == "Watanabe_only"),
    "Intersection of genes with finite MAGMA P values in both preliminary Jansen and Watanabe jobs",
    length(bg_ids),
    "Jansen and Watanabe public analyses substantially overlap in UK Biobank and are not independent replications",
    "Homo sapiens; unranked query; custom background; GO BP/MF/CC, KEGG, Reactome, WikiPathways; g:Profiler g:SCS <0.05; ShinyGO BH FDR <0.05"
  ), stringsAsFactors = FALSE
)
utils::write.table(manifest, file.path(out, "legacy_jansen_watanabe_28_manifest.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
cat("Arquivos históricos prontos em: ", out, "\n", sep = "")
print(manifest)
