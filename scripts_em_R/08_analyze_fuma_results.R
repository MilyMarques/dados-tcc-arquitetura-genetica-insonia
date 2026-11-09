#!/usr/bin/env Rscript
# Análise comparativa reprodutível dos resultados SNP2GENE/MAGMA do FUMA.
# Usa somente R base.

options(stringsAsFactors = FALSE, scipen = 999)

arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(arg)) sub("^--file=", "", arg[1]) else "R/08_analyze_fuma_results.R"
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
input_dir <- file.path(root, "projeto_insonia_tcc", "Resultados")
out_dir <- file.path(input_dir, "analise_padronizada")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

studies <- data.frame(
  study_id = c("hammerschlag_2017", "lane_2019_extreme", "jansen_2019_public", "watanabe_2022_ukb"),
  label = c("Hammerschlag 2017", "Lane 2019: sintomas extremos", "Jansen 2019: UKB público", "Watanabe 2022: UKB público"),
  evidence_unit = c("hammerschlag", "lane", "ukb_public", "ukb_public"),
  zip = c("FUMA_job758083_hammerschlag.zip", "FUMA_job758591_lane.zip", "FUMA_job758172_jansen.zip", "FUMA_job758173_watanabe.zip")
)

read_zip_table <- function(zip_path, member) {
  listed <- utils::unzip(zip_path, list = TRUE)$Name
  if (!member %in% listed) stop("Ausente: ", member, " em ", basename(zip_path))
  td <- tempfile("fuma_"); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  file <- utils::unzip(zip_path, files = member, exdir = td, junkpaths = TRUE)
  utils::read.delim(file[1], check.names = FALSE, na.strings = c("", "NA"), quote = "", stringsAsFactors = FALSE)
}

write_tsv <- function(x, file) utils::write.table(x, file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
safe_num <- function(x) suppressWarnings(as.numeric(x))

summary_rows <- list(); leads <- list(); loci <- list(); mapped <- list(); magma_all <- list()
for (i in seq_len(nrow(studies))) {
  st <- studies[i, ]; zip_path <- file.path(input_dir, st$zip)
  if (!file.exists(zip_path)) stop("ZIP não encontrado: ", zip_path)
  indsig <- read_zip_table(zip_path, "IndSigSNPs.txt")
  lead <- read_zip_table(zip_path, "leadSNPs.txt")
  locus <- read_zip_table(zip_path, "GenomicRiskLoci.txt")
  genes <- read_zip_table(zip_path, "genes.txt")
  eqtl <- read_zip_table(zip_path, "eqtl.txt")
  ci <- read_zip_table(zip_path, "ci.txt")
  magma <- read_zip_table(zip_path, "magma.genes.out")

  magma$P <- safe_num(magma$P)
  tested <- sum(is.finite(magma$P)); bonf <- 0.05 / tested
  magma$bonferroni_significant <- is.finite(magma$P) & magma$P < bonf
  magma$FDR_BH <- NA_real_
  good <- is.finite(magma$P)
  magma$FDR_BH[good] <- p.adjust(magma$P[good], method = "BH")
  magma$FDR_significant <- is.finite(magma$FDR_BH) & magma$FDR_BH < 0.05
  magma$study_id <- st$study_id; magma$label <- st$label; magma$evidence_unit <- st$evidence_unit
  magma_all[[i]] <- magma

  lead$study_id <- st$study_id; lead$label <- st$label; leads[[i]] <- lead
  locus$study_id <- st$study_id; locus$label <- st$label; loci[[i]] <- locus
  genes$study_id <- st$study_id; genes$label <- st$label; genes$evidence_unit <- st$evidence_unit
  # FUMA registra ausência de SNPs posicional/eQTL como 0; somente contagens > 0 são evidência.
  genes$positional_mapping <- !is.na(safe_num(genes$posMapSNPs)) & safe_num(genes$posMapSNPs) > 0
  genes$eqtl_mapping <- !is.na(safe_num(genes$eqtlMapSNPs)) & safe_num(genes$eqtlMapSNPs) > 0
  genes$chromatin_interaction_mapping <- !is.na(genes$ciMap) & genes$ciMap == "Yes"
  mapped[[i]] <- genes

  summary_rows[[i]] <- data.frame(
    study_id = st$study_id, label = st$label, evidence_unit = st$evidence_unit,
    independent_significant_snps = nrow(indsig), lead_snps = nrow(lead), genomic_risk_loci = nrow(locus),
    mapped_genes = nrow(genes), eqtl_pairs = nrow(eqtl), chromatin_interaction_pairs = nrow(ci),
    magma_genes_tested = tested, magma_bonferroni_threshold = bonf,
    magma_bonferroni_significant = sum(magma$bonferroni_significant), magma_fdr_significant = sum(magma$FDR_significant)
  )
}

summary_table <- do.call(rbind, summary_rows)
lead_table <- do.call(rbind, leads)
locus_table <- do.call(rbind, loci)
mapping_table <- do.call(rbind, mapped)
magma_table <- do.call(rbind, magma_all)

# Jansen e Watanabe são outputs UKB sobrepostos: contam como uma unidade na tabela consenso.
magma_sig <- magma_table[magma_table$bonferroni_significant & !is.na(magma_table$SYMBOL) & magma_table$SYMBOL != "", ]
# GENE é o identificador Ensembl único; SYMBOL pode se repetir em anotações históricas.
by_gene <- split(magma_sig, magma_sig$GENE)
consensus <- do.call(rbind, lapply(by_gene, function(d) data.frame(
  GENE = d$GENE[1], SYMBOL = d$SYMBOL[1], n_fuma_studies = length(unique(d$study_id)),
  n_independent_units = length(unique(d$evidence_unit)),
  evidence_units = paste(sort(unique(d$evidence_unit)), collapse = ";"),
  studies = paste(sort(unique(d$study_id)), collapse = ";"), best_magma_p = min(d$P, na.rm = TRUE)
)))
consensus <- consensus[order(-consensus$n_independent_units, -consensus$n_fuma_studies, consensus$best_magma_p), ]

j <- magma_table[magma_table$study_id == "jansen_2019_public", c("GENE", "P")]
w <- magma_table[magma_table$study_id == "watanabe_2022_ukb", c("GENE", "P")]
names(j)[2] <- "P_jansen"; names(w)[2] <- "P_watanabe"
jw <- merge(j, w, by = "GENE")
jw <- jw[is.finite(jw$P_jansen) & is.finite(jw$P_watanabe), ]
jw_summary <- data.frame(comparison = "Jansen UKB vs Watanabe UKB", genes_in_common = nrow(jw),
  spearman_rho_minus_log10_p = cor(-log10(jw$P_jansen), -log10(jw$P_watanabe), method = "spearman"))

write_tsv(summary_table, file.path(out_dir, "01_study_summary.tsv"))
write_tsv(lead_table, file.path(out_dir, "02_lead_snps.tsv"))
write_tsv(locus_table, file.path(out_dir, "03_genomic_risk_loci.tsv"))
write_tsv(mapping_table, file.path(out_dir, "04_mapped_genes_and_evidence.tsv"))
write_tsv(magma_table, file.path(out_dir, "05_magma_all_genes.tsv"))
write_tsv(magma_sig, file.path(out_dir, "06_magma_bonferroni_genes.tsv"))
write_tsv(consensus, file.path(out_dir, "07_magma_consensus_candidates.tsv"))
write_tsv(jw_summary, file.path(out_dir, "08_jansen_watanabe_technical_concordance.tsv"))

cat("Análise concluída em: ", out_dir, "\n", sep = "")
print(summary_table)
print(jw_summary)
