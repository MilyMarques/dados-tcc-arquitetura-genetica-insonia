"""Compara associações curadas do GWAS Catalog para insomnia com loci FUMA principais.

Esta é uma comparação externa descritiva de literatura. Não é uma nova análise GWAS,
não altera o core de 16 genes e não usa associações catalogadas como summary statistics.
"""
from pathlib import Path
import csv
from collections import defaultdict
import zipfile
import io

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "projeto_insonia_tcc"
ANALYSIS = BASE / "Resultados" / "analise_padronizada"


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(rows, path, fields):
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def uniq_join(values):
    return "; ".join(sorted({str(x).strip() for x in values if str(x).strip() and str(x).strip() != "NA"}))


def main():
    catalog = next(BASE.glob("gwas-association-downloaded*"))
    core = read_tsv(ANALYSIS / "11_provisional_consensus_protein_coding.tsv")
    core_symbols = {r["SYMBOL"] for r in core}
    zips = {
        "hammerschlag_2017": "FUMA_job758083_hammerschlag.zip",
        "lane_2019_extreme": "FUMA_job758591_lane.zip",
        "jansen_2019_public": "FUMA_job758172_jansen.zip",
        "watanabe_2022_ukb": "FUMA_job758173_watanabe.zip",
    }
    # Identificadores rsID não dependem de build; usamos SNPs candidatos FUMA (em LD)
    # para atribuir uma associação catalogada a um locus. Não comparamos CHR:POS, pois
    # GWAS Catalog e FUMA podem estar em GRCh38 e GRCh37, respectivamente.
    fuma_snp_to_loci = defaultdict(list)
    for study, name in zips.items():
        with zipfile.ZipFile(BASE / "Resultados" / name) as z:
            with z.open("snps.txt") as fh:
                for r in csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8"), delimiter="\t"):
                    rsid = r.get("rsID", "").upper()
                    if rsid:
                        fuma_snp_to_loci[rsid].append({"study_id": study, "locus_id": r.get("GenomicLocus", "")})

    # Arquivo baixado do Catalog é tab-delimitado, apesar da extensão usualmente .tsv.
    with catalog.open(encoding="utf-8-sig", newline="") as fh:
        assoc = list(csv.DictReader(fh, delimiter="\t"))
    hits = []
    for a in assoc:
        snp_field = (a.get("SNPS", "") + ";" + a.get("STRONGEST SNP-RISK ALLELE", "")).upper()
        rsids = {token.split("-")[0].strip() for token in snp_field.replace(",", ";").split(";") if token.strip().startswith("RS")}
        matching = [l for rsid in rsids for l in fuma_snp_to_loci.get(rsid, [])]
        if not matching:
            continue
        for l in matching:
            reported = a.get("REPORTED GENE(S)", "")
            mapped = a.get("MAPPED_GENE", "")
            gene_tokens = {g.strip() for field in (reported, mapped) for g in field.replace(",", ";").split(";") if g.strip()}
            hits.append({
                "fuma_study_id": l["study_id"], "fuma_locus_id": l["locus_id"], "catalog_rsid_match": "; ".join(sorted(rsids.intersection(fuma_snp_to_loci))),
                "catalog_snp": a.get("SNPS", ""), "catalog_p": a.get("P-VALUE", ""),
                "catalog_study_accession": a.get("STUDY ACCESSION", ""), "catalog_first_author": a.get("FIRST AUTHOR", ""),
                "catalog_study": a.get("STUDY", ""), "reported_trait": a.get("DISEASE/TRAIT", ""),
                "mapped_trait": a.get("MAPPED_TRAIT", ""), "reported_genes": reported, "mapped_genes": mapped,
                "core_gene_in_catalog_annotation": "TRUE" if core_symbols.intersection(gene_tokens) else "FALSE",
                "core_genes_in_catalog_annotation": "; ".join(sorted(core_symbols.intersection(gene_tokens))),
            })

    detail_fields = list(hits[0]) if hits else ["fuma_study_id", "fuma_locus_id", "catalog_rsid_match"]
    write_tsv(hits, ANALYSIS / "32_gwas_catalog_insomnia_x_primary_fuma_loci.tsv", detail_fields)

    grouped = defaultdict(list)
    for h in hits: grouped[(h["fuma_study_id"], h["fuma_locus_id"])].append(h)
    summary = []
    for (study, locus), rows in grouped.items():
        summary.append({
            "fuma_study_id": study, "fuma_locus_id": locus,
            "catalog_associations_in_locus": len(rows),
            "catalog_studies": uniq_join(r["catalog_study_accession"] for r in rows),
            "catalog_authors": uniq_join(r["catalog_first_author"] for r in rows),
            "reported_traits": uniq_join(r["reported_trait"] for r in rows),
            "mapped_traits": uniq_join(r["mapped_trait"] for r in rows),
            "reported_genes": uniq_join(r["reported_genes"] for r in rows),
            "mapped_genes": uniq_join(r["mapped_genes"] for r in rows),
            "core_genes_named_by_catalog": uniq_join(r["core_genes_in_catalog_annotation"] for r in rows),
        })
    summary.sort(key=lambda r: (-int(r["catalog_associations_in_locus"]), r["fuma_study_id"], str(r["fuma_locus_id"])))
    summary_fields = list(summary[0]) if summary else ["fuma_study_id", "fuma_locus_id"]
    write_tsv(summary, ANALYSIS / "33_gwas_catalog_insomnia_locus_summary.tsv", summary_fields)

    metadata = [
        {"item": "analysis_role", "value": "External descriptive concordance with curated GWAS Catalog insomnia associations"},
        {"item": "catalog_file", "value": catalog.name},
        {"item": "catalog_associations_input", "value": len(assoc)},
        {"item": "matching_rule", "value": "Exact rsID overlap between a curated GWAS Catalog association and a FUMA candidate SNP; candidate SNP identifies the FUMA locus. This avoids coordinate/build mismatch."},
        {"item": "important_limit", "value": "GWAS Catalog reported/mapped genes use study-specific curation and are not causal assignments or independent replication"},
    ]
    write_tsv(metadata, ANALYSIS / "34_gwas_catalog_comparison_readme.tsv", ["item", "value"])
    print(f"Catalog associations read: {len(assoc)}; locus-study matches: {len(hits)}; distinct FUMA loci with Catalog evidence: {len(summary)}")


if __name__ == "__main__":
    main()
