#!/usr/bin/env Rscript
# Integra MAGMA competitive gene-set analysis dos dois jobs UKB OpenGWAS.
# Sensibilidade técnica: os datasets compartilham UK Biobank e não redefinem o core.

options(stringsAsFactors = FALSE, scipen = 999)
arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- if (length(arg)) sub("^--file=", "", arg[1]) else "R/19_integrate_opengwas_ukb_gene_sets.R"
root <- normalizePath(file.path(dirname(script), ".."), winslash="/", mustWork=TRUE)
out_dir <- file.path(root,"projeto_insonia_tcc","Resultados","analise_padronizada")
zip_dir <- file.path(root,"projeto_insonia_tcc","Resultados_opengwas")
studies <- data.frame(study_id=c("ukb_a_13","ukb_b_3957"), label=c("OpenGWAS ukb-a-13 (Neale Lab)","OpenGWAS ukb-b-3957 (PHESANT)"), zip=c("FUMA_job761601_UKB-a-13.zip","FUMA_job761609_UKB-b-3957.zip"), stringsAsFactors=FALSE)

read_gsa <- function(zip_path) {
  td <- tempfile("gsa_"); dir.create(td); on.exit(unlink(td,recursive=TRUE),add=TRUE)
  members <- utils::unzip(zip_path,list=TRUE)$Name; member <- members[members=="magma.gsa.out"]
  if(!length(member)) stop("magma.gsa.out ausente: ",basename(zip_path))
  f <- utils::unzip(zip_path,files=member[1],exdir=td,junkpaths=TRUE)
  utils::read.table(f,header=TRUE,comment.char="#",fill=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
}
write_tsv <- function(x,n) utils::write.table(x,file.path(out_dir,n),sep="\t",row.names=FALSE,quote=FALSE,na="")

all <- list(); sums <- list()
for(i in seq_len(nrow(studies))) {
  s <- studies[i,]; x <- read_gsa(file.path(zip_dir,s$zip)); x$P <- suppressWarnings(as.numeric(x$P)); valid <- is.finite(x$P)
  x$FDR_BH <- NA_real_; x$FDR_BH[valid] <- p.adjust(x$P[valid],"BH")
  x$bonferroni_significant <- valid & x$P < 0.05/sum(valid); x$FDR_significant <- valid & x$FDR_BH < 0.05
  x$study_id <- s$study_id; x$label <- s$label; all[[i]] <- x
  sums[[i]] <- data.frame(study_id=s$study_id,label=s$label,gene_sets_tested=sum(valid),bonferroni_threshold=0.05/sum(valid),n_bonferroni=sum(x$bonferroni_significant),n_fdr=sum(x$FDR_significant),best_p=min(x$P[valid]),role="UKB sensitivity only",stringsAsFactors=FALSE)
}
all <- do.call(rbind,all); sums <- do.call(rbind,sums)
by_set <- split(all,all$FULL_NAME)
compare <- do.call(rbind,lapply(by_set,function(d) data.frame(FULL_NAME=d$FULL_NAME[1],TYPE=d$TYPE[1],NGENES=median(d$NGENES),n_ukb_runs_fdr=sum(d$FDR_significant),n_ukb_runs_bonferroni=sum(d$bonferroni_significant),best_p=min(d$P),best_fdr=min(d$FDR_BH),runs_fdr=paste(d$study_id[d$FDR_significant],collapse=";"),stringsAsFactors=FALSE)))
compare <- compare[order(-compare$n_ukb_runs_fdr,compare$best_fdr),]
reportable <- compare[compare$n_ukb_runs_fdr>=1,]; reportable$interpretation <- ifelse(reportable$n_ukb_runs_fdr==2,"Reproduzido nas duas sensibilidades UKB (não independente)","Apenas uma sensibilidade UKB")
write_tsv(sums,"35_opengwas_ukb_magma_gene_set_study_summary.tsv")
write_tsv(all,"36_opengwas_ukb_magma_gene_set_all.tsv")
write_tsv(compare,"37_opengwas_ukb_magma_gene_set_compare.tsv")
write_tsv(reportable,"38_opengwas_ukb_magma_gene_set_fdr_reportable.tsv")
write_tsv(data.frame(item=c("role","method","limit"),value=c("Sensitivity only","MAGMA competitive gene-set; BH FDR within each UKB job","Both datasets are UK Biobank based and not independent validation")),"39_opengwas_ukb_magma_gene_set_readme.tsv")
cat("UKB OpenGWAS gene-set analysis concluded. FDR reportable sets: ",nrow(reportable),"\n",sep="")
