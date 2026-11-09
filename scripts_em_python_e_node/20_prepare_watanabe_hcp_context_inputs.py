"""Prepara lista contextual de Watanabe 2022 para enriquecimento/STRING secundário.

Entrada: Tabela suplementar 14 (3.898 genes FUMA candidatos) e Tabela 37
(genes de alta confiança). HCP=1 identifica os 289 genes priorizados pelos autores.
Não altera o core16 nem é validação independente.
"""
from pathlib import Path
import csv
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "projeto_insonia_tcc"
SUPP = BASE / "DADOS SUPLEMENTARES - Watanabe 2022"
OUT = BASE / "enrichment_inputs" / "watanabe_2022_hcp289_context"
OUT.mkdir(parents=True, exist_ok=True)


def rows(ws, header_row=2):
    head = [str(x).strip() if x is not None else "" for x in next(ws.iter_rows(min_row=header_row, max_row=header_row, values_only=True))]
    return [dict(zip(head, row)) for row in ws.iter_rows(min_row=header_row + 1, values_only=True) if row and row[0] is not None]


def write_lines(values, name):
    (OUT / name).write_text("\n".join(str(x) for x in values if x not in (None, "")) + "\n", encoding="utf-8")


def write_headered(values, name, header="GENE"):
    # O GENE2FUNC começa a ler na segunda linha de arquivos enviados.
    (OUT / name).write_text(header + "\n" + "\n".join(str(x) for x in values if x not in (None, "")) + "\n", encoding="utf-8")


def main():
    xlsx = next(SUPP.glob("Supplementary Tables*.xlsx"))
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    candidates = rows(wb["Supplementary Table 14"])
    hc = rows(wb["Supplementary Table 37"])
    # A tabela 14 tem 3.898 registros gene-locus; depois de remover repetições de
    # ENSG, há 3.526 genes únicos.
    bg = {r["ENSG"]: r for r in candidates}
    # HCP=1 é a prioridade publicada: 289 genes, conforme artigo/tabela suplementar.
    hcp = {r["ENSG ID"]: r for r in hc if r.get("HCP") == 1}
    if len(hcp) != 289:
        raise RuntimeError(f"Esperados 289 HCP genes; encontrados {len(hcp)}")
    missing = set(hcp) - set(bg)
    # Um HCP não está na lista de genes ingenuamente mapeados. Para que a
    # própria lista de interesse esteja contida no universo testável, o
    # background contextual é a união (3.899 genes), com essa exceção
    # explicitamente registrada no manifesto.
    bg.update({k: {"ENSG": k, "Symbol": hcp[k].get("Symbol", "")} for k in missing})

    hcp_rows = [hcp[k] for k in sorted(hcp)]
    bg_rows = [bg[k] for k in sorted(bg)]
    write_lines([r["ENSG ID"] for r in hcp_rows], "watanabe_hcp289_ensembl_ids.txt")
    write_lines([r["Symbol"] for r in hcp_rows], "watanabe_hcp289_gene_symbols.txt")
    write_lines([r["ENSG"] for r in bg_rows], "watanabe_fuma_mapped_unique_background_ensembl_ids.txt")
    write_lines([r["Symbol"] for r in bg_rows], "watanabe_fuma_mapped_unique_background_gene_symbols.txt")
    write_headered([r["ENSG ID"] for r in hcp_rows], "gene2func_watanabe_hcp289_ensembl.tsv")
    write_headered([r["ENSG"] for r in bg_rows], "gene2func_watanabe_hcp289_background_ensembl.tsv")
    with (OUT / "watanabe_hcp289_gene_table.tsv").open("w", newline="", encoding="utf-8") as fh:
        fields = ["ENSG ID", "Symbol", "HUGO", "chr", "Locus", "Label", "HCP"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(hcp_rows)
    manifest = [
        ("analysis_role", "External literature-context analysis based on Watanabe 2022; not a new discovery and does not redefine the core16"),
        ("input_definition", "289 HCP genes (HCP=1) from Supplementary Table 37"),
        ("input_gene_count", str(len(hcp))),
        ("background_definition", "3,898 gene-locus records naively mapped by FUMA from the 554 Watanabe risk loci (Supplementary Table 14), representing 3,526 unique ENSG genes; plus 1 HCP gene absent from that table so the query is contained in the background"),
        ("background_gene_count", str(len(bg))),
        ("recommended_enrichment", "g:Profiler/ShinyGO with HCP289 vs the unique FUMA-mapped gene background; GO BP/MF/CC, KEGG, Reactome, WikiPathways"),
        ("hcp_not_in_table14", ";".join(sorted(missing)) if missing else "none"),
        ("important_limit", "Watanabe includes UK Biobank/23andMe meta-analysis; the public UKB component overlaps the current project. Interpret as published literature context, not independent validation."),
    ]
    with (OUT / "watanabe_hcp289_manifest.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["item", "value"]); w.writerows(manifest)
    print(f"HCP genes: {len(hcp)}; mapped-gene background: {len(bg)}; output: {OUT}")


if __name__ == "__main__":
    main()
