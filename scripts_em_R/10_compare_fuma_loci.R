#!/usr/bin/env Rscript
# Compara formalmente loci FUMA por interseção de intervalos GRCh37.

options(stringsAsFactors = FALSE, scipen = 999)
arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(arg)) sub("^--file=", "", arg[1]) else "R/10_compare_fuma_loci.R"
root <- normalizePath(file.path(dirname(script), ".."), winslash = "/", mustWork = TRUE)
dir <- file.path(root, "projeto_insonia_tcc", "Resultados", "analise_padronizada")
loci <- utils::read.delim(file.path(dir, "03_genomic_risk_loci.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, name) utils::write.table(x, file.path(dir, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "")

loci$chr <- as.integer(loci$chr); loci$start <- as.integer(loci$start); loci$end <- as.integer(loci$end)
loci$locus_key <- paste(loci$study_id, loci$GenomicLocus, sep = "__")
loci$evidence_unit <- ifelse(loci$study_id %in% c("jansen_2019_public", "watanabe_2022_ukb"), "ukb_public", loci$study_id)

# Cada par só é chamado de sobreposto quando seus intervalos FUMA no mesmo cromossomo se intersectam >= 1 bp.
edges <- list(); k <- 1
for (i in seq_len(nrow(loci) - 1)) for (j in (i + 1):nrow(loci)) {
  if (loci$study_id[i] == loci$study_id[j] || loci$chr[i] != loci$chr[j]) next
  overlap_start <- max(loci$start[i], loci$start[j]); overlap_end <- min(loci$end[i], loci$end[j])
  if (overlap_start <= overlap_end) {
    overlap_bp <- overlap_end - overlap_start + 1
    edges[[k]] <- data.frame(
      study_a = loci$study_id[i], locus_a = loci$GenomicLocus[i], lead_a = loci$rsID[i], p_a = loci$p[i], start_a = loci$start[i], end_a = loci$end[i],
      study_b = loci$study_id[j], locus_b = loci$GenomicLocus[j], lead_b = loci$rsID[j], p_b = loci$p[j], start_b = loci$start[j], end_b = loci$end[j],
      chr = loci$chr[i], overlap_start = overlap_start, overlap_end = overlap_end, overlap_bp = overlap_bp,
      reciprocal_overlap_a = overlap_bp / (loci$end[i] - loci$start[i] + 1),
      reciprocal_overlap_b = overlap_bp / (loci$end[j] - loci$start[j] + 1), stringsAsFactors = FALSE
    ); k <- k + 1
  }
}
edges <- if (length(edges)) do.call(rbind, edges) else data.frame()

# Componentes conexos: um cluster pode ligar >= 2 estudos mesmo que nem todo par se sobreponha diretamente.
parent <- seq_len(nrow(loci))
find_root <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
join <- function(a, b) { ra <- find_root(a); rb <- find_root(b); if (ra != rb) parent[rb] <<- ra }
if (nrow(edges)) for (r in seq_len(nrow(edges))) {
  a <- match(paste(edges$study_a[r], edges$locus_a[r], sep = "__"), loci$locus_key)
  b <- match(paste(edges$study_b[r], edges$locus_b[r], sep = "__"), loci$locus_key)
  join(a, b)
}
clusters <- vapply(seq_len(nrow(loci)), find_root, integer(1))
loci$locus_cluster <- paste0("CL", match(clusters, unique(clusters)))

cluster_rows <- lapply(split(loci, loci$locus_cluster), function(d) data.frame(
  locus_cluster = d$locus_cluster[1], chr = d$chr[1], cluster_start = min(d$start), cluster_end = max(d$end),
  n_fuma_loci = nrow(d), n_studies = length(unique(d$study_id)), n_independent_units = length(unique(d$evidence_unit)),
  studies = paste(sort(unique(d$study_id)), collapse = ";"), evidence_units = paste(sort(unique(d$evidence_unit)), collapse = ";"),
  lead_rsids = paste(d$rsID, collapse = ";"), best_p = min(as.numeric(d$p)), stringsAsFactors = FALSE
))
cluster_table <- do.call(rbind, cluster_rows)
cluster_table <- cluster_table[order(-cluster_table$n_independent_units, cluster_table$best_p), ]

if (nrow(edges)) {
  pair_summary <- aggregate(overlap_bp ~ study_a + study_b, edges, length)
  names(pair_summary)[3] <- "overlapping_locus_pairs"
} else pair_summary <- data.frame()

write_tsv(edges, "13_pairwise_fuma_locus_overlaps.tsv")
write_tsv(loci, "14_locus_cluster_members.tsv")
write_tsv(cluster_table, "15_cross_study_locus_clusters.tsv")
write_tsv(pair_summary, "16_pairwise_locus_overlap_summary.tsv")
cat("Comparação de loci concluída em: ", dir, "\n", sep = "")
print(pair_summary)
