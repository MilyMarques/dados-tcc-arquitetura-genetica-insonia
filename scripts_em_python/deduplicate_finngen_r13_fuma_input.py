"""Remove duplicatas globais CHR:BP:A1:A2 apos liftover FinnGen R13->hg19.

Mantem a primeira ocorrencia no arquivo ordenado em GRCh38; gera auditoria JSON.
"""
from __future__ import annotations
import gzip, json, sqlite3, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "projeto_insonia_tcc" / "fuma_inputs_external" / "finngen_r13_f5_insomnia"
PATH = OUTDIR / "finngen_r13_f5_insomnia_full.FUMA.txt.gz"
PARTIAL = OUTDIR / ".finngen_r13_f5_insomnia_full.dedup.partial.txt.gz"

def main() -> None:
    tmp = tempfile.NamedTemporaryFile(suffix=".sqlite", delete=False); db_path=Path(tmp.name); tmp.close()
    db=sqlite3.connect(db_path); db.execute("PRAGMA journal_mode=OFF"); db.execute("PRAGMA synchronous=OFF"); db.execute("CREATE TABLE seen(k TEXT PRIMARY KEY)")
    report={"input":str(PATH),"rows_input":0,"rows_kept":0,"duplicates_removed":0,"rule":"Mantida a primeira ocorrencia de CHR:BP:A1:A2 apos liftover GRCh38->hg19."}
    try:
        with gzip.open(PATH,"rt",encoding="utf-8") as src, gzip.open(PARTIAL,"wt",encoding="utf-8",newline="") as dst:
            dst.write(src.readline())
            for line in src:
                report["rows_input"] += 1
                x=line.rstrip("\n").split("\t"); key=f"{x[1]}:{x[2]}:{x[3]}:{x[4]}"
                try: db.execute("INSERT INTO seen VALUES (?)",(key,))
                except sqlite3.IntegrityError: report["duplicates_removed"] += 1; continue
                dst.write(line); report["rows_kept"] += 1
        db.commit(); PARTIAL.replace(PATH)
    finally:
        db.close(); db_path.unlink(missing_ok=True); PARTIAL.unlink(missing_ok=True)
    (OUTDIR/"finngen_r13_f5_insomnia.post_liftover_dedup.json").write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(report,indent=2))
if __name__=="__main__": main()
