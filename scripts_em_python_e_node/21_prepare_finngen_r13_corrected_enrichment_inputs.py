"""Cria inputs de enrichment diretamente do ZIP FUMA FinnGen R13 corrigido.

O script preserva os inputs da primeira rodada e escreve uma nova pasta para o
job 763972 (FinnGen R13 corrigido/full). A query é composta por genes
protein-coding significativos no MAGMA após Bonferroni; o background são todos
os genes protein-coding com P MAGMA finito desse mesmo job.
"""
from pathlib import Path
import csv
import io
import math
import zipfile


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "projeto_insonia_tcc"
ZIP = BASE / "Resultados_finn" / "FUMA_job763972_finnr13_corrigido-full.zip"
OUT = BASE / "enrichment_inputs" / "finngen_r13_corrigido_external"


def read_tsv_from_zip(zf, member):
    with zf.open(member) as fh:
        return list(csv.DictReader(io.TextIOWrapper(fh, encoding="utf-8"), delimiter="\t"))


def write_tsv(path, values):
    path.write_text("GENE\n" + "\n".join(values) + "\n", encoding="utf-8")


def write_txt(path, values):
    path.write_text("\n".join(values) + "\n", encoding="utf-8")


def main():
    if not ZIP.exists():
        raise FileNotFoundError(f"ZIP FUMA corrigido não encontrado: {ZIP}")
    OUT.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(ZIP) as zf:
        magma = read_tsv_from_zip(zf, "magma.genes.out")

    # O arquivo magma.genes.out já é o universo de genes testados pelo MAGMA.
    # Não usar genes.txt aqui: ele contém apenas genes mapeados pelos loci SNP2GENE,
    # não a lista genome-wide testada pelo MAGMA.
    records = []
    for r in magma:
        try:
            p = float(r["P"])
        except (KeyError, ValueError, TypeError):
            continue
        if not math.isfinite(p):
            continue
        symbol = r.get("SYMBOL", "").strip()
        if not symbol or symbol == "NA":
            continue
        records.append({"ensg": r["GENE"], "symbol": symbol, "p": p})

    records.sort(key=lambda x: x["ensg"])
    threshold = 0.05 / len(records)
    query = [r for r in records if r["p"] < threshold]

    write_tsv(OUT / "finngen_r13_corrigido_magma_bonferroni_ensembl.tsv", [r["ensg"] for r in query])
    write_tsv(OUT / "finngen_r13_corrigido_magma_bonferroni_symbols.tsv", [r["symbol"] for r in query])
    write_tsv(OUT / "finngen_r13_corrigido_magma_testable_background_ensembl.tsv", [r["ensg"] for r in records])
    write_tsv(OUT / "finngen_r13_corrigido_magma_testable_background_symbols.tsv", [r["symbol"] for r in records])
    # Versões sem cabeçalho: apropriadas para colar nos campos web.
    write_txt(OUT / "finngen_r13_corrigido_query_ensembl_paste.txt", [r["ensg"] for r in query])
    write_txt(OUT / "finngen_r13_corrigido_query_symbols_paste.txt", [r["symbol"] for r in query])
    write_txt(OUT / "finngen_r13_corrigido_background_ensembl_paste.txt", [r["ensg"] for r in records])
    write_txt(OUT / "finngen_r13_corrigido_background_symbols_paste.txt", [r["symbol"] for r in records])

    with (OUT / "finngen_r13_corrigido_enrichment_manifest.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["item", "value"])
        w.writerows([
            ("source_job", "FUMA_job763972_finnr13_corrigido-full.zip"),
            ("analysis_role", "External Finnish clinical-insomnia validation; does not redefine the four-study core16"),
            ("query_definition", "Genes with symbol and MAGMA P < 0.05 / number of MAGMA-testable genes with symbol"),
            ("query_gene_count", len(query)),
            ("background_definition", "Genes with finite MAGMA P and non-empty symbol in the corrected FinnGen R13 FUMA job"),
            ("background_gene_count", len(records)),
            ("bonferroni_threshold", f"{threshold:.16g}"),
            ("recommended_platforms", "GENE2FUNC and g:Profiler use Ensembl; ShinyGO uses symbols; always provide the custom background"),
        ])
    print(f"Criados: query={len(query)}, background={len(records)}, threshold={threshold:.6g}")
    print(OUT)


if __name__ == "__main__":
    main()
