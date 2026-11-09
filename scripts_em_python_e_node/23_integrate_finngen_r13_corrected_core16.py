"""Atualiza a validação FinnGen do core16 usando exclusivamente o job FUMA corrigido."""
from pathlib import Path
import csv
import io
import math
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "projeto_insonia_tcc"
OUT = BASE / "Resultados" / "analise_padronizada"
ZIP = BASE / "Resultados_finn" / "FUMA_job763972_finnr13_corrigido-full.zip"


def read_zip_tsv(member):
    with zipfile.ZipFile(ZIP) as zf, zf.open(member) as fh:
        return list(csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8"), delimiter="\t"))


def write(rows, name):
    with (OUT / name).open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "22_core16_final_evidence_table.tsv").open(encoding="utf-8") as fh:
        core = list(csv.DictReader(fh, delimiter="\t"))
    magma = read_zip_tsv("magma.genes.out")
    by_gene = {r["GENE"]: r for r in magma}
    finite = [float(r["P"]) for r in magma if r.get("P") and math.isfinite(float(r["P"]))]
    threshold = 0.05 / len(finite)
    rows = []
    for c in core:
        r = by_gene.get(c["GENE"])
        p = float(r["P"]) if r else math.nan
        rows.append({
            "GENE": c["GENE"], "SYMBOL": c["SYMBOL"], "priority_tier": c["priority_tier"],
            "finngen_r13_corrected_magma_P": "" if not math.isfinite(p) else f"{p:.12g}",
            "finngen_r13_corrected_bonferroni_threshold": f"{threshold:.12g}",
            "finngen_r13_corrected_bonferroni_significant": str(math.isfinite(p) and p < threshold).upper(),
            "interpretation": "Suporte externo FinnGen R13 corrigido" if math.isfinite(p) and p < threshold else "Sem suporte Bonferroni FinnGen R13 corrigido",
        })
    write(rows, "49_core16_x_finngen_r13_corrigido_magma.tsv")
    summary = [{"metric": "source_job", "value": ZIP.name},
               {"metric": "magma_genes_tested", "value": str(len(finite))},
               {"metric": "bonferroni_threshold", "value": f"{threshold:.12g}"},
               {"metric": "core16_tested", "value": str(len(rows))},
               {"metric": "core16_bonferroni_supported", "value": str(sum(x["finngen_r13_corrected_bonferroni_significant"] == "TRUE" for x in rows))}]
    write(summary, "50_core16_finngen_r13_corrigido_validation_summary.tsv")
    print(f"FinnGen corrigido: {len(finite)} genes MAGMA; limiar={threshold:.6g}; core16 Bonf={summary[-1]['value']}")


if __name__ == "__main__":
    main()
