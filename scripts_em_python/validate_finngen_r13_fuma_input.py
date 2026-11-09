"""Valida o input FinnGen R13 preparado para FUMA, inclusive duplicatas globais.

Uso:
  py projeto_insonia_tcc/pipeline_fuma/validate_finngen_r13_fuma_input.py --full
"""
from __future__ import annotations
import argparse, gzip, json, re, sqlite3, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PATH = ROOT / "projeto_insonia_tcc" / "fuma_inputs_external" / "finngen_r13_f5_insomnia" / "finngen_r13_f5_insomnia_full.FUMA.txt.gz"
RSID, DNA = re.compile(r"^rs\d+$", re.I), {"A", "C", "G", "T"}

def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("--full", action="store_true"); args = ap.parse_args()
    r = {"study_id":"finngen_r13_f5_insomnia", "file":str(PATH), "pass":False, "errors":[], "rows":0, "p_min":None, "gws_variants":0, "duplicates_global":"not_checked"}
    db = db_path = None
    if args.full:
        tmp = tempfile.NamedTemporaryFile(suffix=".sqlite", delete=False); db_path = Path(tmp.name); tmp.close()
        db = sqlite3.connect(db_path); db.execute("PRAGMA journal_mode=OFF"); db.execute("PRAGMA synchronous=OFF"); db.execute("CREATE TABLE seen(k TEXT PRIMARY KEY)")
    try:
        with gzip.open(PATH, "rt", encoding="utf-8") as f:
            expected = ["SNP","CHR","BP","A1","A2","BETA","SE","P","N"]
            if f.readline().rstrip("\n").split("\t") != expected: r["errors"].append("cabecalho FUMA inesperado")
            for line_no, line in enumerate(f, 2):
                x = line.rstrip("\n").split("\t")
                try:
                    snp, chrom, bp, a1, a2, beta, se, p, n = x; chrom, bp, n = int(chrom), int(bp), int(n); beta, se, p = float(beta), float(se), float(p)
                    valid = len(x)==9 and 1<=chrom<=22 and bp>0 and n==497916 and RSID.fullmatch(snp) and a1 in DNA and a2 in DNA and a1!=a2 and se>0 and 0<p<=1
                except (ValueError, TypeError): valid=False
                if not valid: r["errors"].append(f"linha invalida: {line_no}"); break
                r["rows"] += 1; r["p_min"] = p if r["p_min"] is None else min(r["p_min"],p); r["gws_variants"] += p < 5e-8
                if db:
                    try: db.execute("INSERT INTO seen VALUES (?)", (f"{chrom}:{bp}:{a1}:{a2}",))
                    except sqlite3.IntegrityError: r["errors"].append("duplicata global detectada"); break
        r["duplicates_global"] = "checked" if args.full else "not_checked"; r["pass"] = bool(r["rows"] and not r["errors"])
    finally:
        if db: db.close()
        if db_path: db_path.unlink(missing_ok=True)
    out=PATH.with_suffix("").with_suffix(".validation.json"); out.write_text(json.dumps(r,indent=2,ensure_ascii=False),encoding="utf-8")
    print("PASS" if r["pass"] else "FAIL", json.dumps(r,indent=2)); raise SystemExit(not r["pass"])
if __name__ == "__main__": main()
