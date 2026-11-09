"""Prepara a camada global do GWAS Catalog para análise descritiva.

Usa somente SNP_GENE_IDS (IDs Ensembl fornecidos pelo próprio GWAS Catalog),
evitando tratar a coluna textual MAPPED_GENE como identificação inequívoca.
Os arquivos resultantes são apropriados para g:Profiler, GENE2FUNC e STRING.
Não fazem atribuição causal de genes.
"""
from pathlib import Path
import csv
import re

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "projeto_insonia_tcc"
SOURCE = BASE / "gwas-association-downloaded_2026-08-12-MONDO_0013600_insomnia.tsv"
OUT = BASE / "enrichment_inputs" / "gwas_catalog_insomnia_global"

ENS_RE = re.compile(r"^ENSG\d{11}$")


def split(value):
    if not value or value.strip().upper() in {"NR", "NA"}:
        return []
    return [x.strip() for x in value.split(";") if x.strip()]


def uniq(values):
    return sorted(set(v for v in values if v))


def write_lines(name, values, header=None):
    text = ((header + "\n") if header else "") + "\n".join(values) + "\n"
    (OUT / name).write_text(text, encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with SOURCE.open(encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    by_gene = {}
    raw_labels = []
    for r in rows:
        ids = [x for x in split(r.get("SNP_GENE_IDS", "")) if ENS_RE.fullmatch(x)]
        labels = split(r.get("MAPPED_GENE", ""))
        raw_labels.extend(labels)
        for ens in ids:
            d = by_gene.setdefault(ens, {"ENSG": ens, "n_catalog_associations": 0, "rsids": [], "study_accessions": [],
                                         "authors": [], "reported_traits": [], "mapped_traits": [], "mapped_gene_labels": [], "p_values": []})
            d["n_catalog_associations"] += 1
            d["rsids"].extend(split(r.get("SNPS", "")))
            d["study_accessions"].extend(split(r.get("STUDY ACCESSION", "")))
            d["authors"].extend(split(r.get("FIRST AUTHOR", "")))
            d["reported_traits"].extend(split(r.get("DISEASE/TRAIT", "")))
            d["mapped_traits"].extend(split(r.get("MAPPED_TRAIT", "")))
            d["mapped_gene_labels"].extend(labels)
            try:
                d["p_values"].append(float(r.get("P-VALUE", "")))
            except ValueError:
                pass

    agg = []
    for ens, d in sorted(by_gene.items()):
        agg.append({"ENSG": ens, "n_catalog_associations": d["n_catalog_associations"],
                    "n_distinct_rsids": len(uniq(d["rsids"])), "n_distinct_studies": len(uniq(d["study_accessions"])),
                    "min_catalog_p": min(d["p_values"]) if d["p_values"] else "",
                    "rsids": "; ".join(uniq(d["rsids"])), "study_accessions": "; ".join(uniq(d["study_accessions"])),
                    "authors": "; ".join(uniq(d["authors"])), "reported_traits": "; ".join(uniq(d["reported_traits"])),
                    "mapped_traits": "; ".join(uniq(d["mapped_traits"])),
                    "mapped_gene_labels_in_same_associations": "; ".join(uniq(d["mapped_gene_labels"]))})

    ensg = [r["ENSG"] for r in agg]
    write_lines("gwas_catalog_insomnia_global_ensembl_ids.txt", ensg)
    write_lines("gene2func_gwas_catalog_insomnia_global_ensembl.tsv", ensg, header="GENE")
    write_lines("gwas_catalog_insomnia_global_string_ensembl_ids.txt", ensg)
    write_lines("gwas_catalog_insomnia_raw_mapped_gene_labels.txt", uniq(raw_labels))
    with (OUT / "gwas_catalog_insomnia_global_gene_provenance.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(agg[0]), delimiter="\t"); w.writeheader(); w.writerows(agg)
    manifest = [
        ("analysis_role", "Descriptive global GWAS Catalog layer; not causal gene prioritisation and not a replacement for core16"),
        ("source_file", SOURCE.name),
        ("catalog_association_rows", len(rows)),
        ("catalog_studies", len(uniq([r.get("STUDY ACCESSION", "") for r in rows]))),
        ("global_gene_definition", "Unique valid Ensembl IDs in SNP_GENE_IDS, as provided by GWAS Catalog association records"),
        ("global_gene_count", len(ensg)),
        ("raw_mapped_gene_label_count", len(uniq(raw_labels))),
        ("recommended_gprofiler", "Query Ensembl IDs; default annotated human background; Benjamini-Hochberg FDR < 0.01 to mirror Lopes 2021; All results off"),
        ("recommended_string", "Use Ensembl query; organism Homo sapiens; no added interactors; high confidence >=0.7; retain mapping report and show connected component separately if desired"),
        ("critical_limit", "Catalog associations are curated top associations and include study/phenotype/locus repetition. Mapped genes are annotations, not causal assignments."),
    ]
    with (OUT / "gwas_catalog_insomnia_global_manifest.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["item", "value"]); w.writerows(manifest)
    print(f"Catalog rows={len(rows)}; unique ENSG={len(ensg)}; raw gene labels={len(uniq(raw_labels))}")


if __name__ == "__main__":
    main()
