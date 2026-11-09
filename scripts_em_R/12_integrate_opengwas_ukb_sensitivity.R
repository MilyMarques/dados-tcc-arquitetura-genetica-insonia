#!/usr/bin/env Rscript
# Integra os dois SNP2GENE/MAGMA OpenGWAS-UKB como análise de sensibilidade.
# Não altera o consenso principal, pois as amostras são UK Biobank sobrepostas.

options(stringsAsFactors = FALSE, scipen = 999)
arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(arg)) sub("^--file=", "", arg[1]) else "R/12_integrate_opengwas_ukb_sensitivity.R"
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
main_dir <- file.path(root, "projeto_insonia_tcc", "Resultados", "analise_padronizada")
opengwas_dir <- file.path(root, "projeto_insonia_tcc", "Resultados_opengwas")

studies <- data.frame(
  study_id = c("ukb_a_13", "ukb_b_3957"),
  label = c("OpenGWAS ukb-a-13 (Neale Lab)", "OpenGWAS ukb-b-3957 (PHESANT)"),
  zip = c("FUMA_job761601_UKB-a-13.zip", "FUMA_job761609_UKB-b-3957.zip"),
  stringsAsFactors = FALSE
)

read_zip_table <- function(zip_path, member) {
  listed <- utils::unzip(zip_path, list = TRUE)$Name
  if (!member %in% listed) stop("Ausente: ", member, " em ", basename(zip_path))
  td <- tempfile("fuma_"); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  file <- utils::unzip(zip_path, files = member, exdir = td, junkpaths = TRUE)
  utils::read.delim(file[1], check.names = FALSE, na.strings = c("", "NA"), quote = "", stringsAsFactors = FALSE)
}
write_tsv <- function(x, name) utils::write.table(x, file.path(main_dir, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "")

core <- utils::read.delim(file.path(main_dir, "11_provisional_consensus_protein_coding.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
core <- core[, c("GENE", "SYMBOL", "priority_tier", "n_independent_units_magma_bonf", "n_mapping_modalities", "best_magma_p")]

all_magma <- list(); summary_rows <- list()
for (i in seq_len(nrow(studies))) {
  st <- studies[i, ]; zip_path <- file.path(opengwas_dir, st$zip)
  if (!file.exists(zip_path)) stop("ZIP não encontrado: ", zip_path)
  magma <- read_zip_table(zip_path, "magma.genes.out")
  magma$P <- suppressWarnings(as.numeric(magma$P))
  tested <- sum(is.finite(magma$P)); bonf <- 0.05 / tested
  magma$FDR_BH <- NA_real_; valid <- is.finite(magma$P)
  magma$FDR_BH[valid] <- p.adjust(magma$P[valid], method = "BH")
  magma$bonferroni_significant <- valid & magma$P < bonf
  magma$FDR_significant <- is.finite(magma$FDR_BH) & magma$FDR_BH < 0.05
  magma$study_id <- st$study_id; magma$label <- st$label
  all_magma[[i]] <- magma
  summary_rows[[i]] <- data.frame(study_id = st$study_id, label = st$label, magma_genes_tested = tested,
    magma_bonferroni_threshold = bonf, magma_bonferroni_significant = sum(magma$bonferroni_significant),
    magma_fdr_significant = sum(magma$FDR_significant), role = "UKB sensitivity only", stringsAsFactors = FALSE)
}

magma_all <- do.call(rbind, all_magma)
study_summary <- do.call(rbind, summary_rows)

# Uma linha por gene do núcleo e por estudo UKB, incluindo ausência de teste se ela ocorresse.
grid <- merge(core, studies[, c("study_id", "label")], all = TRUE)
ukb_core <- merge(grid, magma_all[, c("GENE", "SYMBOL", "study_id", "P", "FDR_BH", "bonferroni_significant", "FDR_significant", "NSNPS", "N")], by = c("GENE", "study_id"), all.x = TRUE, suffixes = c("_core", "_ukb"))
ukb_core$SYMBOL <- ifelse(!is.na(ukb_core$SYMBOL_core), ukb_core$SYMBOL_core, ukb_core$SYMBOL_ukb)
ukb_core$SYMBOL_core <- NULL; ukb_core$SYMBOL_ukb <- NULL
ukb_core$validated_magma_bonferroni <- !is.na(ukb_core$bonferroni_significant) & ukb_core$bonferroni_significant
ukb_core$validated_magma_fdr <- !is.na(ukb_core$FDR_significant) & ukb_core$FDR_significant
ukb_core <- ukb_core[order(ukb_core$priority_tier, ukb_core$SYMBOL, ukb_core$study_id), ]

by_gene <- split(ukb_core, ukb_core$GENE)
ukb_summary <- do.call(rbind, lapply(by_gene, function(d) data.frame(
  GENE = d$GENE[1], SYMBOL = d$SYMBOL[1], priority_tier = d$priority_tier[1],
  n_ukb_sensitivity_runs_tested = sum(!is.na(d$P)),
  n_ukb_sensitivity_runs_bonferroni = sum(d$validated_magma_bonferroni),
  n_ukb_sensitivity_runs_fdr = sum(d$validated_magma_fdr),
  ukb_runs_bonferroni = paste(sort(d$study_id[d$validated_magma_bonferroni]), collapse = ";"),
  best_ukb_magma_p = if (any(is.finite(d$P))) min(d$P, na.rm = TRUE) else NA_real_,
  sensitivity_interpretation = "Sensibilidade técnica UKB; não adiciona unidade de evidência independente.",
  stringsAsFactors = FALSE
)))
ukb_summary <- ukb_summary[order(-ukb_summary$n_ukb_sensitivity_runs_bonferroni, ukb_summary$best_ukb_magma_p), ]

write_tsv(study_summary, "19_opengwas_ukb_magma_summary.tsv")
write_tsv(ukb_core, "20_core16_x_opengwas_ukb_magma.tsv")
write_tsv(ukb_summary, "21_core16_opengwas_ukb_sensitivity_summary.tsv")
cat("Integração UKB OpenGWAS concluída em: ", main_dir, "\n", sep = "")
print(study_summary)
print(ukb_summary[, c("SYMBOL", "n_ukb_sensitivity_runs_bonferroni", "best_ukb_magma_p")])
