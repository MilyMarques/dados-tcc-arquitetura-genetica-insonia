"""Prepara summary statistics de insônia para FUMA SNP2GENE/MAGMA.

Este script lê arquivos grandes em blocos. Arquivos brutos nunca são alterados.
Os outputs são candidatos a upload: devem passar pelo validador companion antes
de serem enviados ao FUMA.

Uso (a partir da pasta projeto_codex):
    python projeto_insonia_tcc/pipeline_fuma/process_sumstats_to_fuma.py --study all
    python projeto_insonia_tcc/pipeline_fuma/process_sumstats_to_fuma.py --study watanabe_2022_ukb
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "projeto_insonia_tcc" / "fuma_inputs_final"
CHUNK_SIZE = 500_000
INFO_MIN = 0.80
MAF_MIN = 0.001
ALLOWED_CHR = set(range(1, 23))  # autossomos: protocolo único para comparação
P_GWS = 5e-8


@dataclass(frozen=True)
class StudySpec:
    key: str
    input_path: Path
    compression: str | None
    sep: str
    columns: dict[str, str]
    effect_type: str  # "beta" ou "or"
    fixed_n: int | None
    include_snp: bool
    require_rsid: bool
    note: str


SPECS = {
    "hammerschlag_2017": StudySpec(
        key="hammerschlag_2017",
        input_path=PROJECT_ROOT / "Hammerschlag_2017" /
        "Hammerschlag_NatGenet2017_insomnia_sumstats-full_090617.txt.gz",
        compression="gzip",
        sep="\t",
        columns={"snp": "SNP", "chr": "CHR", "bp": "POS", "a1": "A1", "a2": "A2",
                 "beta": "BETA", "se": "SE", "p": "P", "n": "N_analyzed",
                 "eaf": "EAF", "info": "INFO"},
        effect_type="beta", fixed_n=None, include_snp=True, require_rsid=True,
        note="Arquivo original GRCh37; BETA original é usado porque reproduz P com SE.",
    ),
    "jansen_2019_public": StudySpec(
        key="jansen_2019_public",
        input_path=PROJECT_ROOT / "JansenPR_2019" / "Insomnia_sumstats_Jansenetal.txt.gz",
        compression="gzip",
        sep="\t",
        columns={"snp": "SNP", "chr": "CHR", "bp": "BP", "a1": "A1", "a2": "A2",
                 "or": "OR", "se": "SE", "p": "P", "n": "N", "maf": "MAF", "info": "INFO"},
        effect_type="or", fixed_n=None, include_snp=True, require_rsid=True,
        note="Versão pública UK Biobank; BETA é calculado como log(OR).",
    ),
    "lane_2019_extreme": StudySpec(
        key="lane_2019_extreme",
        input_path=PROJECT_ROOT / "Lane_2019" / "Lane_30804566_fullUKBB_Insomnia_summary_stats.zip",
        compression="zip",
        sep=r"\s+",
        columns={"chr": "CHR", "bp": "BP", "a1": "ALLELE1", "a2": "ALLELE0",
                 "beta": "BETA_INSOMNIA", "se": "SE_INSOMNIA", "p": "P_INSOMNIA",
                 "eaf": "A1FREQ", "info": "INFO"},
        effect_type="beta", fixed_n=237627, include_snp=False, require_rsid=False,
        note="Insomnia extreme symptoms versus controls. Sem coluna SNP para evitar IDs chr:pos não-rsID; FUMA usa CHR+BP.",
    ),
    "watanabe_2022_ukb": StudySpec(
        key="watanabe_2022_ukb",
        input_path=PROJECT_ROOT / "Watanabe_2022" /
        "insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz",
        compression="gzip",
        sep="\t",
        columns={"snp": "RSID_UKB", "chr": "CHR", "bp": "BP", "a1": "A1", "a2": "A2",
                 "or": "OR", "se": "SE", "p": "P", "n": "NMISS", "maf": "MAF",
                 "info": "INFO_UKB"},
        effect_type="or", fixed_n=None, include_snp=True, require_rsid=True,
        note="Componente UK Biobank. RSID_UKB é usado como SNP; A1/A2 são os alelos do modelo e BETA=log(OR).",
    ),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def clean_chr(series: pd.Series) -> pd.Series:
    value = series.astype("string").str.upper().str.strip().str.replace("^CHR", "", regex=True)
    value = value.replace({"X": "23", "Y": "24", "MT": "25", "M": "25"})
    return pd.to_numeric(value, errors="coerce")


def add_count(report: dict[str, int], name: str, before: int, after: int) -> None:
    report[name] += before - after


def valid_dna_allele(series: pd.Series) -> pd.Series:
    # Mantém SNVs e indels representados apenas por A/C/G/T; FUMA aceita indels
    # se os dois alelos forem informados exatamente.
    return series.str.fullmatch(r"[ACGT]+", na=False)


def standardize_chunk(raw: pd.DataFrame, spec: StudySpec, report: dict[str, int]) -> pd.DataFrame:
    missing = sorted(set(spec.columns.values()) - set(raw.columns))
    if missing:
        raise ValueError(f"{spec.key}: colunas ausentes: {missing}")

    use = raw.loc[:, list(spec.columns.values())].rename(
        columns={source: target.upper() for target, source in spec.columns.items()}
    ).copy()
    report["raw_rows"] += len(use)

    # Textos
    for col in ("A1", "A2", "SNP"):
        if col in use:
            use[col] = use[col].astype("string").str.upper().str.strip()

    # Números
    use["CHR"] = clean_chr(use["CHR"])
    for col in ("BP", "BETA", "OR", "SE", "P", "N", "EAF", "MAF", "INFO"):
        if col in use:
            use[col] = pd.to_numeric(use[col], errors="coerce")
    if "N" not in use:
        use["N"] = spec.fixed_n

    required = ["CHR", "BP", "A1", "A2", "SE", "P", "N"]
    required.append("BETA" if spec.effect_type == "beta" else "OR")
    if spec.include_snp:
        required.append("SNP")
    before = len(use)
    use = use.dropna(subset=required)
    add_count(report, "missing_required", before, len(use))

    before = len(use)
    use = use.loc[(use["P"] > 0) & (use["P"] <= 1) & (use["SE"] > 0) & (use["N"] > 0)]
    if spec.effect_type == "or":
        use = use.loc[use["OR"] > 0]
        use["BETA"] = np.log(use["OR"])
    else:
        use = use.loc[np.isfinite(use["BETA"])]
    add_count(report, "invalid_statistics", before, len(use))

    before = len(use)
    use = use.loc[use["CHR"].isin(ALLOWED_CHR) & (use["BP"] > 0)]
    add_count(report, "non_autosomal_or_invalid_position", before, len(use))

    before = len(use)
    use = use.loc[valid_dna_allele(use["A1"]) & valid_dna_allele(use["A2"]) &
                  (use["A1"] != use["A2"])]
    add_count(report, "invalid_alleles", before, len(use))

    if "INFO" in use and use["INFO"].notna().any():
        before = len(use)
        use = use.loc[use["INFO"] >= INFO_MIN]
        add_count(report, "info_below_0_8", before, len(use))

    if "MAF" not in use and "EAF" in use:
        use["MAF"] = np.minimum(use["EAF"], 1 - use["EAF"])
    if "MAF" in use and use["MAF"].notna().any():
        before = len(use)
        use = use.loc[(use["MAF"] >= MAF_MIN) & (use["MAF"] <= 0.5)]
        add_count(report, "maf_below_0_001_or_invalid", before, len(use))

    if spec.include_snp and spec.require_rsid:
        before = len(use)
        use = use.loc[use["SNP"].str.fullmatch(r"RS\d+", case=False, na=False)]
        add_count(report, "missing_or_non_rsid", before, len(use))

    use["CHR"] = use["CHR"].astype("int32")
    use["BP"] = use["BP"].astype("int64")
    use["N"] = use["N"].astype("int64")
    before_dedup = len(use)
    use = use.sort_values("P", kind="stable").drop_duplicates(
        subset=["CHR", "BP", "A1", "A2"], keep="first"
    )
    report["chunk_duplicate_variants_removed"] += before_dedup - len(use)
    report["rows_after_qc"] += len(use)
    return use


def process_study(spec: StudySpec) -> Path:
    if not spec.input_path.exists():
        raise FileNotFoundError(spec.input_path)
    out_dir = OUTPUT_ROOT / spec.key
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / f"{spec.key}.FUMA.txt.gz"
    partial = out_dir / f".{spec.key}.partial.txt.gz"
    qc_path = out_dir / f"{spec.key}.qc.json"
    schema_path = out_dir / f"{spec.key}.submission_metadata.json"
    partial.unlink(missing_ok=True)

    report = {
        "study_id": spec.key, "input_file": str(spec.input_path), "input_bytes": spec.input_path.stat().st_size,
        "build": "GRCh37", "allowed_chromosomes": "1-22", "info_min": INFO_MIN,
        "maf_min": MAF_MIN, "chunksize": CHUNK_SIZE, "raw_rows": 0, "rows_after_qc": 0,
        "missing_required": 0, "invalid_statistics": 0, "non_autosomal_or_invalid_position": 0,
        "invalid_alleles": 0, "info_below_0_8": 0, "maf_below_0_001_or_invalid": 0,
        "missing_or_non_rsid": 0, "chunk_duplicate_variants_removed": 0,
        "started_at": utc_now(), "note": spec.note,
    }
    header = True
    output_columns = (["SNP"] if spec.include_snp else []) + ["CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"]

    try:
        with gzip.open(partial, "wt", encoding="utf-8", newline="") as handle:
            reader = pd.read_csv(
                spec.input_path, sep=spec.sep, compression=spec.compression, dtype="string",
                chunksize=CHUNK_SIZE, low_memory=False,
            )
            for number, raw in enumerate(reader, start=1):
                cleaned = standardize_chunk(raw, spec, report)
                cleaned.loc[:, output_columns].to_csv(handle, sep="\t", index=False, header=header)
                header = False
                print(f"[{spec.key}] bloco {number}: {len(raw):,} lidas; {len(cleaned):,} mantidas", flush=True)
        if header:
            raise RuntimeError(f"{spec.key}: nenhuma variante passou pelo QC.")
        partial.replace(output)
    except Exception:
        partial.unlink(missing_ok=True)
        raise

    report["finished_at"] = utc_now()
    report["fuma_candidate_file"] = str(output)
    report["fuma_candidate_bytes"] = output.stat().st_size
    report["output_columns"] = output_columns
    qc_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    schema_path.write_text(json.dumps({
        "study_id": spec.key, "fuma_file": output.name, "build": "GRCh37",
        "columns": output_columns, "effect_allele": "A1", "non_effect_allele": "A2",
        "p_column": "P", "sample_size_column": "N", "candidate_status": "VALIDATE_BEFORE_UPLOAD",
        "note": spec.note,
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", default="all", choices=["all", *SPECS.keys()])
    args = parser.parse_args()
    selected = SPECS.values() if args.study == "all" else [SPECS[args.study]]
    for spec in selected:
        print(f"\n=== PROCESSANDO {spec.key} ===", flush=True)
        file = process_study(spec)
        print(f"Candidato FUMA gerado: {file}", flush=True)
    print("\nPróximo passo obrigatório: execute validate_fuma_inputs.py --full antes do upload.")


if __name__ == "__main__":
    main()
