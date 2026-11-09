"""Validação independente dos inputs FUMA produzidos a partir de VCF OpenGWAS."""
from __future__ import annotations

import argparse, gzip, json, re, sqlite3, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INPUT_ROOT = ROOT / "projeto_insonia_tcc" / "fuma_inputs_opengwas"
STUDIES = ("finn_b_f5_insomnia", "ukb_b_3957", "ukb_a_13")
RSID = re.compile(r"^rs\d+$", re.IGNORECASE)
DNA = {"A", "C", "G", "T"}


def validate(study: str, full: bool) -> dict:
    path = INPUT_ROOT / study / f"{study}.FUMA.txt.gz"
    result = {"study_id": study, "file": str(path), "pass": False, "errors": [], "rows": 0,
              "p_min": None, "gws_variants": 0, "duplicates_global": "not_checked"}
    db = None; db_path = None
    if full:
        tmp = tempfile.NamedTemporaryFile(prefix=f"{study}_", suffix=".sqlite", delete=False)
        db_path = Path(tmp.name); tmp.close(); db = sqlite3.connect(db_path)
        db.execute("PRAGMA journal_mode=OFF"); db.execute("PRAGMA synchronous=OFF"); db.execute("CREATE TABLE seen(k TEXT PRIMARY KEY)")
    try:
        with gzip.open(path, "rt", encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            expected = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"]
            if header != expected:
                result["errors"].append(f"cabeçalho inesperado: {header}")
            for line_number, line in enumerate(f, start=2):
                values = line.rstrip("\n").split("\t")
                if len(values) != 9:
                    result["errors"].append(f"linha {line_number}: número de colunas inválido"); break
                snp, chrom, bp, a1, a2, beta, se, p, n = values
                try:
                    chrom_i, bp_i, n_i = int(chrom), int(bp), int(n)
                    beta_f, se_f, p_f = float(beta), float(se), float(p)
                except ValueError:
                    result["errors"].append(f"linha {line_number}: valor numérico inválido"); break
                valid = chrom_i in range(1, 23) and bp_i > 0 and n_i > 0 and se_f > 0 and 0 < p_f <= 1 and RSID.fullmatch(snp) and a1 in DNA and a2 in DNA and a1 != a2
                if not valid:
                    result["errors"].append(f"linha {line_number}: valores/alelos inválidos"); break
                result["rows"] += 1; result["p_min"] = min(result["p_min"] if result["p_min"] is not None else 1.0, p_f)
                result["gws_variants"] += int(p_f < 5e-8)
                if db is not None:
                    key = f"{chrom}:{bp}:{a1}:{a2}"
                    try: db.execute("INSERT INTO seen(k) VALUES (?)", (key,))
                    except sqlite3.IntegrityError:
                        result["errors"].append("duplicata global detectada"); break
        result["duplicates_global"] = "checked" if full else "not_checked"
        result["pass"] = bool(result["rows"] > 0 and not result["errors"])
    finally:
        if db is not None: db.close()
        if db_path is not None: db_path.unlink(missing_ok=True)
    out = path.with_suffix("").with_suffix(".validation.json")
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", choices=("all", *STUDIES), default="all")
    parser.add_argument("--full", action="store_true", help="verifica duplicatas globais com SQLite")
    args = parser.parse_args(); selected = STUDIES if args.study == "all" else (args.study,)
    results = [validate(study, args.full) for study in selected]
    for r in results: print(r["study_id"], "PASS" if r["pass"] else "FAIL", r)
    if not all(r["pass"] for r in results): raise SystemExit(1)


if __name__ == "__main__":
    main()
