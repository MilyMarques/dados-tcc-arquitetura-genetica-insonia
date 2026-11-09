#!/usr/bin/env Rscript
# Integra MAGMA competitive gene-set analyses dos quatro FUMA SNP2GENE principais.
# Esta análise usa o sinal gene-based de todo o genoma; não é enrichment da lista core16.

options(stringsAsFactors = FALSE, scipen = 999)
arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(arg)) sub("^--file=", "", arg[1]) else "R/16_integrate_magma_gene_sets.R"
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "projeto_insonia_tcc", "Resultados", "analise_padronizada")
res_dir <- file.path(root, "projeto_insonia_tcc", "Resultados")

studies <- data.frame(
  study_id = c("hammerschlag_2017", "lane_2019_extreme", "jansen_2019_public", "watanabe_2022_ukb"),
  label = c("Hammerschlag 2017", "Lane 2019", "Jansen 2019 public UKB", "Watanabe 2022 public UKB"),
  zip = c("FUMA_job758083_hammerschlag.zip", "FUMA_job758591_lane.zip", "FUMA_job758172_jansen.zip", "FUMA_job758173_watanabe.zip"),
  unit = c("hammerschlag", "lane", "ukb_public", "ukb_public"), stringsAsFactors = FALSE
)

read_gsa <- function(zip_path) {
  members <- utils::unzip(zip_path, list = TRUE)$Name
  member <- members[members == "magma.gsa.out"]
  if (!length(member)) stop("magma.gsa.out ausente em ", basename(zip_path))
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  utils::unzip(zip_path, files = member[1], exdir = dirname(tmp), junkpaths = TRUE)
  extracted <- file.path(dirname(tmp), basename(member[1]))
  x <- utils::read.table(extracted, header = TRUE, comment.char = "#", fill = TRUE,
                         stringsAsFactors = FALSE, check.names = FALSE)
  unlink(extracted)
  if (!("P" %in% names(x)) || !("FULL_NAME" %in% names(x))) stop("Formato MAGMA GSA inesperado: ", basename(zip_path))
  x$P <- suppressWarnings(as.numeric(x$P)); x
}
write_tsv <- function(x, name) utils::write.table(x, file.path(out_dir, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "")

all <- list(); sums <- list()
for (i in seq_len(nrow(studies))) {
  s <- studies[i, ]; x <- read_gsa(file.path(res_dir, s$zip))
  valid <- is.finite(x$P)
  x$FDR_BH <- NA_real_; x$FDR_BH[valid] <- p.adjust(x$P[valid], method = "BH")
  x$bonf_threshold <- 0.05 / sum(valid)
  x$bonferroni_significant <- valid & x$P < x$bonf_threshold
  x$FDR_significant <- valid & x$FDR_BH < 0.05
  x$study_id <- s$study_id; x$label <- s$label; x$independent_unit <- s$unit
  all[[i]] <- x
  sums[[i]] <- data.frame(study_id=s$study_id,label=s$label,independent_unit=s$unit,
    gene_sets_tested=sum(valid),bonferroni_threshold=0.05/sum(valid),
    n_bonferroni=sum(x$bonferroni_significant),n_fdr=sum(x$FDR_significant),
    best_p=min(x$P[valid]),stringsAsFactors=FALSE)
}
all <- do.call(rbind, all); study_summary <- do.call(rbind, sums)

# Uma linha por gene-set: consistência entre estudos. Jansen e Watanabe contam
# como uma unidade UKB, usando o melhor p dentro da unidade para a descrição.
sets <- split(all, all$FULL_NAME)
consensus <- do.call(rbind, lapply(sets, function(d) {
  unit_rows <- split(d, d$independent_unit)
  best_by_unit <- lapply(unit_rows, function(z) z[which.min(z$P), , drop=FALSE])
  b <- do.call(rbind, best_by_unit)
  data.frame(
    FULL_NAME=d$FULL_NAME[1], TYPE=d$TYPE[1], NGENES=median(d$NGENES, na.rm=TRUE),
    n_studies_tested=nrow(d), n_independent_units_tested=nrow(b),
    n_studies_fdr=sum(d$FDR_significant), n_independent_units_fdr=sum(b$FDR_significant),
    n_studies_bonferroni=sum(d$bonferroni_significant), n_independent_units_bonferroni=sum(b$bonferroni_significant),
    best_p=min(d$P,na.rm=TRUE), best_fdr=min(d$FDR_BH,na.rm=TRUE),
    units_fdr=paste(b$independent_unit[b$FDR_significant],collapse=";"),
    studies_fdr=paste(d$study_id[d$FDR_significant],collapse=";"),
    stringsAsFactors=FALSE)
}))
consensus <- consensus[order(-consensus$n_independent_units_fdr, consensus$best_fdr, consensus$best_p), ]

# Conjunto de resultados interpretáveis: FDR em >=1 estudo; marca robustez se >=2 unidades.
reportable <- consensus[consensus$n_independent_units_fdr >= 1, ]
reportable$interpretation <- ifelse(reportable$n_independent_units_fdr >= 2,
  "Consistente em >=2 unidades independentes", "Significativo em uma unidade; exploratório")

write_tsv(study_summary, "25_magma_gene_set_study_summary.tsv")
write_tsv(all, "26_magma_gene_set_all_studies.tsv")
write_tsv(consensus, "27_magma_gene_set_cross_study_summary.tsv")
write_tsv(reportable, "28_magma_gene_set_fdr_reportable.tsv")

readme <- data.frame(
  item=c("analysis", "method", "multiple_testing", "independence_rule", "reporting_rule", "scope"),
  value=c(
    "MAGMA competitive gene-set analysis integrated across four primary FUMA jobs",
    "Tests whether genes in a set show stronger gene-based GWAS association than other genes, controlling MAGMA covariates",
    "Benjamini-Hochberg FDR and Bonferroni recalculated within each study over all finite-P gene sets",
    "Jansen and Watanabe treated as one UKB public unit for cross-study recurrence",
    "FDR in >=1 unit is reportable; >=2 independent units is labeled consistent",
    "Genome-wide gene-set analysis; distinct from over-representation analyses of the core16 gene list"
  ),stringsAsFactors=FALSE)
write_tsv(readme, "29_magma_gene_set_readme.tsv")
cat("MAGMA gene-set integration complete: ", out_dir, "\n", sep="")
print(study_summary)
cat("Reportable FDR gene sets: ",nrow(reportable),"\n",sep="")
