"""Valida os candidatos FUMA gerados por process_sumstats_to_fuma.py.

Uso:
    python projeto_insonia_tcc/pipeline_fuma/validate_fuma_inputs.py --study all --full

--full também verifica duplicatas globais por CHR:BP:A1:A2 em um banco SQLite
temporário no disco. É mais lento, mas é a validação final antes do upload.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sqlite3
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
INPUT_ROOT = ROOT / "projeto_insonia_tcc" / "fuma_inputs_final"
STUDIES = ("hammerschlag_2017", "jansen_2019_public", "lane_2019_extreme", "watanabe_2022_ukb")
REQUIRED = {"CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"}


def validate(study: str, full_duplicates: bool) -> dict:
    path = INPUT_ROOT / study / f"{study}.FUMA.txt.gz"
    if not path.exists():
        raise FileNotFoundError(path)
    summary = {"study_id": study, "file": str(path), "pass": False, "errors": [],
               "rows": 0, "p_min": None, "gws_variants": 0, "duplicates_global": None}
    seen = None
    db_path = None
    if full_duplicates:
        tmp = tempfile.NamedTemporaryFile(prefix=f"{study}_", suffix=".sqlite", delete=False)
        db_path = Path(tmp.name); tmp.close()
        seen = sqlite3.connect(db_path)
        seen.execute("PRAGMA journal_mode=OFF")
        seen.execute("PRAGMA synchronous=OFF")
        seen.execute("CREATE TABLE seen (k TEXT PRIMARY KEY)")
    try:
        for chunk in pd.read_csv(path, sep="\t", compression="gzip", chunksize=500_000):
            missing = REQUIRED - set(chunk.columns)
            if missing:
                summary["errors"].append(f"colunas ausentes: {sorted(missing)}")
                break
            chr_ok = pd.to_numeric(chunk["CHR"], errors="coerce").isin(range(1, 23))
            bp = pd.to_numeric(chunk["BP"], errors="coerce")
            p = pd.to_numeric(chunk["P"], errors="coerce")
            n = pd.to_numeric(chunk["N"], errors="coerce")
            allele = chunk["A1"].astype("string").str.fullmatch(r"[ACGT]+", na=False) & \
                     chunk["A2"].astype("string").str.fullmatch(r"[ACGT]+", na=False) & \
                     (chunk["A1"] != chunk["A2"])
            ok = chr_ok & (bp > 0) & (p > 0) & (p <= 1) & (n > 0) & allele
            if not bool(ok.all()):
                summary["errors"].append(f"{int((~ok).sum())} linhas inválidas")
            summary["rows"] += len(chunk)
            summary["p_min"] = min(summary["p_min"] if summary["p_min"] is not None else 1.0, float(p.min()))
            summary["gws_variants"] += int((p < 5e-8).sum())
            if seen is not None:
                key = (chunk["CHR"].astype(str) + ":" + chunk["BP"].astype(str) + ":" +
                       chunk["A1"].astype(str) + ":" + chunk["A2"].astype(str))
                before = seen.total_changes
                seen.executemany("INSERT OR IGNORE INTO seen(k) VALUES (?)", ((x,) for x in key))
                if seen.total_changes - before != len(chunk):
                    summary["errors"].append("duplicatas globais detectadas")
        summary["duplicates_global"] = "checked" if full_duplicates else "not_checked"
        summary["pass"] = bool(summary["rows"] > 0 and not summary["errors"])
    finally:
        if seen is not None:
            seen.close()
        if db_path is not None:
            db_path.unlink(missing_ok=True)
    report = path.with_suffix("").with_suffix(".validation.json")
    report.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", choices=("all", *STUDIES), default="all")
    parser.add_argument("--full", action="store_true", help="verifica duplicatas globais")
    args = parser.parse_args()
    selected = STUDIES if args.study == "all" else (args.study,)
    results = [validate(study, args.full) for study in selected]
    for result in results:
        print(result["study_id"], "PASS" if result["pass"] else "FAIL", result)
    if not all(x["pass"] for x in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
