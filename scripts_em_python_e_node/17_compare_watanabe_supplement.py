"""Compara o core de 16 genes com as Tabelas Suplementares 14 e 37 de Watanabe 2022.

Uso (na pasta projeto_codex):
    python scripts/17_compare_watanabe_supplement.py
"""
from pathlib import Path
import csv
import openpyxl
from datetime import timedelta

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "projeto_insonia_tcc"
OUT = BASE / "Resultados" / "analise_padronizada"
SUPP = BASE / "DADOS SUPLEMENTARES - Watanabe 2022"


def read_lines(path):
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def sheet_rows(ws, header_row):
    header = [str(x).strip() if x is not None else "" for x in next(ws.iter_rows(min_row=header_row, max_row=header_row, values_only=True))]
    rows = []
    for row in ws.iter_rows(min_row=header_row + 1, values_only=True):
        if not row or row[0] is None:
            continue
        rows.append(dict(zip(header, row)))
    return rows


def yes(value):
    return "TRUE" if bool(value) else "FALSE"


def excel_scalar(value):
    """Evita que células de locus numeradas pelo Excel sejam exportadas como duração."""
    if isinstance(value, timedelta):
        # O suplemento formata algumas células de GenomicLocus como hora/data;
        # a conversão não é confiável. A presença no FUMA é preservada, mas o
        # identificador daquele locus não é inferido artificialmente.
        return ""
    return value


def main():
    core_file = BASE / "enrichment_inputs" / "core16_ensembl_ids.txt"
    core = read_lines(core_file)
    core_table = OUT / "11_provisional_consensus_protein_coding.tsv"
    with core_table.open(encoding="utf-8", newline="") as fh:
        core_info = {r["GENE"]: r for r in csv.DictReader(fh, delimiter="\t")}

    xlsx = next(SUPP.glob("Supplementary Tables*.xlsx"))
    wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    fuma_rows = sheet_rows(wb["Supplementary Table 14"], 2)
    hc_rows = sheet_rows(wb["Supplementary Table 37"], 2)
    fuma = {r["ENSG"]: r for r in fuma_rows}
    hc = {r["ENSG ID"]: r for r in hc_rows}

    fields = [
        "GENE", "SYMBOL", "priority_tier", "current_magma_units", "current_mapping_modalities",
        "watanabe_fuma_mapped", "watanabe_fuma_locus", "watanabe_fuma_in_locus",
        "watanabe_high_confidence", "watanabe_hc_label", "watanabe_hcp",
        "comparison_interpretation",
    ]
    out_rows = []
    for gene in core:
        c = core_info[gene]
        f = fuma.get(gene, {})
        h = hc.get(gene, {})
        fmap = gene in fuma
        hcp = h.get("HCP") == 1
        hlabel = h.get("Label", "") if h else ""
        if hcp:
            interp = "Concordante: gene core atual e HCP no estudo Watanabe"
        elif fmap:
            interp = "Concordante: gene core atual tambem mapeado por FUMA em Watanabe"
        else:
            interp = "Gene core atual adicional: nao listado no mapeamento FUMA publicado de Watanabe"
        out_rows.append({
            "GENE": gene, "SYMBOL": c["SYMBOL"], "priority_tier": c["priority_tier"],
            "current_magma_units": c["magma_bonf_units"], "current_mapping_modalities": c["n_mapping_modalities"],
            "watanabe_fuma_mapped": yes(fmap), "watanabe_fuma_locus": excel_scalar(f.get("GenomicLocus", "")) if f else "",
            "watanabe_fuma_in_locus": f.get("inLocus", "") if f else "",
            "watanabe_high_confidence": yes(bool(h)), "watanabe_hc_label": hlabel,
            "watanabe_hcp": yes(hcp), "comparison_interpretation": interp,
        })

    output = OUT / "30_core16_x_watanabe_supplement_comparison.tsv"
    with output.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(out_rows)
    summary = OUT / "31_core16_x_watanabe_supplement_summary.tsv"
    with summary.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["core_genes", len(out_rows)])
        writer.writerow(["mapped_by_published_watanabe_fuma", sum(r["watanabe_fuma_mapped"] == "TRUE" for r in out_rows)])
        writer.writerow(["listed_as_watanabe_high_confidence", sum(r["watanabe_high_confidence"] == "TRUE" for r in out_rows)])
        writer.writerow(["watanabe_HCP_genes", sum(r["watanabe_hcp"] == "TRUE" for r in out_rows)])
        writer.writerow(["important_limit", "Watanabe public UKB data overlap with the current analysis; this is literature concordance, not independent replication."])
    print(f"Wrote: {output}\nWrote: {summary}")


if __name__ == "__main__":
    main()
