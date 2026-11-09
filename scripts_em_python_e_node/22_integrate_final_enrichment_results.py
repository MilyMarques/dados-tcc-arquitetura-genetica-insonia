"""Integra as rodadas finais de GENE2FUNC, ShinyGO e g:Profiler.

Produz arquivos em Resultados/analise_padronizada sem alterar arquivos brutos.
O registro separa resultados confirmatórios (p ajustado < 0,05) de saídas
exploratórias, em particular as exportações "All results" do g:Profiler.
"""
from pathlib import Path
import csv
import io
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "projeto_insonia_tcc"
OUT = BASE / "Resultados" / "analise_padronizada"
OUT.mkdir(parents=True, exist_ok=True)

G2F = {
    "core16": (BASE / "Resultados_GENE2FUNC" / "FUMA_gene2func761816_16_genes.zip", "Núcleo principal; quatro estudos", 16),
    "legacy28": (BASE / "Resultados_GENE2FUNC" / "FUMA_gene2func762928_28_genes_legacy.zip", "Sensibilidade histórica Jansen–Watanabe; não independente", 28),
    "watanabe_hcp289": (BASE / "Resultados_GENE2FUNC" / "FUMA_gene2func763975_Watanabe_2022_HCP289_contextual.zip", "Contexto publicado de Watanabe; não redefine core16", 289),
    "finngen_r13_corrected": (BASE / "Resultados_GENE2FUNC" / "FUMA_gene2func764514_finngen_R13.zip", "Validação externa FinnGen R13 corrigido", 21),
}


def g2f_rows(zip_path):
    with zipfile.ZipFile(zip_path) as zf:
        text = io.TextIOWrapper(zf.open("GS.txt"), encoding="utf-8")
        return list(csv.DictReader(text, delimiter="\t"))


def sig_g2f(zip_path):
    return [r for r in g2f_rows(zip_path) if float(r["adjP"]) < 0.05]


def theme(term):
    t = term.lower()
    if any(x in t for x in ("developmental", "hedgehog", "l1cam", "midbrain", "gaba", "neurotype")):
        return "desenvolvimento e organização neuronal"
    if "gap junction" in t or "gap_junction" in t:
        return "comunicação intercelular por junções comunicantes"
    if any(x in t for x in ("ras", "mapk", "growth hormone")):
        return "sinalização intracelular"
    if any(x in t for x in ("autophagy", "post-translational", "fanconi")):
        return "homeostase e manutenção celular"
    if any(x in t for x in ("alzheimer", "parkinson", "huntington", "prion", "amyotrophic", "neurodegeneration")):
        return "conjunto de doença neurológica; interpretar como genes/processos compartilhados"
    return "outro"


def main():
    registry = []
    terms = []
    for analysis_id, (zip_path, role, query_size) in G2F.items():
        sig = sig_g2f(zip_path)
        cats = "; ".join(sorted({r["Category"] for r in sig})) or "nenhuma"
        registry.append({
            "analysis_id": analysis_id, "platform": "FUMA GENE2FUNC", "role": role,
            "query_size": query_size, "custom_background": "sim", "correction": "Benjamini-Hochberg FDR < 0.05",
            "n_significant_gene_sets": len(sig), "significant_categories": cats,
            "status": "significativo" if sig else "sem termos significativos", "source": str(zip_path.relative_to(BASE)),
        })
        for r in sig:
            # GWAS Catalog, TF targets e conjuntos posicionais são anotacionais; a tabela de vias foca categorias funcionais.
            if r["Category"] in {"Canonical_Pathways", "Curated_gene_sets", "KEGG", "Reactome", "Wikipathways", "GO_bp", "Cell_type_signature"}:
                terms.append({"analysis_id": analysis_id, "platform": "FUMA GENE2FUNC", "database_category": r["Category"],
                              "term": r["GeneSet"], "n_overlap": r["N_overlap"], "adjusted_p": r["adjP"],
                              "theme": theme(r["GeneSet"]), "driver_genes": r["genes"],
                              "interpretation_limit": "enrichment de conjunto gênico; não implica causalidade individual"})

    shiny = BASE / "Resultados_ShinyGO" / "ShinyGO_HCP_289_Watanabe" / "enrichment_all.csv"
    with shiny.open(encoding="utf-8-sig", newline="") as fh:
        sh_rows = list(csv.DictReader(fh))
    sh_sig = [r for r in sh_rows if float(r["Enrichment FDR"]) < 0.05]
    registry.append({"analysis_id": "watanabe_hcp289", "platform": "ShinyGO", "role": "Contexto publicado de Watanabe; não redefine core16",
                     "query_size": 289, "custom_background": "sim; background contextual Watanabe informado", "correction": "FDR < 0.05",
                     "n_significant_gene_sets": len(sh_sig), "significant_categories": "KEGG", "status": "significativo",
                     "source": str(shiny.relative_to(BASE))})
    for r in sh_sig:
        terms.append({"analysis_id": "watanabe_hcp289", "platform": "ShinyGO", "database_category": "KEGG",
                      "term": r["Pathway"], "n_overlap": r["nGenes"], "adjusted_p": r["Enrichment FDR"],
                      "theme": theme(r["Pathway"]), "driver_genes": r["Genes"],
                      "interpretation_limit": "enrichment de conjunto gênico; não implica causalidade individual"})

    for analysis_id, role, query_size in [("core16", "Núcleo principal; quatro estudos", 16), ("legacy28", "Sensibilidade histórica Jansen–Watanabe", 28), ("watanabe_hcp289", "Contexto publicado Watanabe", 289), ("finngen_r13_corrected", "Validação externa FinnGen R13 corrigido", 21)]:
        registry.append({"analysis_id": analysis_id, "platform": "g:Profiler g:GOSt", "role": role, "query_size": query_size,
                         "custom_background": "sim", "correction": "g:SCS < 0.05", "n_significant_gene_sets": 0,
                         "significant_categories": "nenhuma", "status": "sem termos significativos; All results é exploratório", "source": "Resultados_gProfiler"})
    # A exportação final do ShinyGO FinnGen não contém FDR < 0,05.
    registry.append({"analysis_id": "finngen_r13_corrected", "platform": "ShinyGO", "role": "Validação externa FinnGen R13 corrigido", "query_size": 21,
                     "custom_background": "sim", "correction": "FDR < 0.05", "n_significant_gene_sets": 0,
                     "significant_categories": "nenhuma", "status": "sem termos significativos", "source": "Resultados_ShinyGO/ShinyGO_finngen_r13/enrichment_all.csv"})
    for analysis_id, role, query_size in [("core16", "Núcleo principal; quatro estudos", 16), ("legacy28", "Sensibilidade histórica Jansen–Watanabe", 28)]:
        registry.append({"analysis_id": analysis_id, "platform": "ShinyGO", "role": role, "query_size": query_size,
                         "custom_background": "sim", "correction": "FDR < 0.05", "n_significant_gene_sets": 0,
                         "significant_categories": "nenhuma", "status": "sem termos significativos", "source": "Resultados_ShinyGO"})

    fields = list(registry[0])
    with (OUT / "47_registro_integrado_analises_funcionais.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t"); w.writeheader(); w.writerows(registry)
    fields = list(terms[0])
    with (OUT / "48_termos_funcionais_significativos_integrados.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t"); w.writeheader(); w.writerows(terms)
    print(f"Registro: {len(registry)} análises; termos significativos registrados: {len(terms)}")


if __name__ == "__main__":
    main()
