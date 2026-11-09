"""Prepara FinnGen R13 F5_INSOMNIA (GRCh38) para FUMA SNP2GENE (hg19).

O arquivo nativo FinnGen e lido em streaming. Antes da exportacao, as
coordenadas GRCh38 sao convertidas para hg19 com a cadeia UCSC local
hg38ToHg19.over.chain.gz. O arquivo final e .txt.gz, com cabecalho FUMA.

Uso no PowerShell, a partir de projeto_codex:
  py projeto_insonia_tcc/pipeline_fuma/process_finngen_r13_to_fuma.py
"""
from __future__ import annotations

import bisect
import gzip
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "projeto_insonia_tcc"
SOURCE = BASE / "sumstats_brutas" / "finngen_r13" / "finngen_R13_F5_INSOMNIA.gz"
CHAIN = BASE / "reference_data" / "hg38ToHg19.over.chain.gz"
OUTDIR = BASE / "fuma_inputs_external" / "finngen_r13_f5_insomnia"
OUT = OUTDIR / "finngen_r13_f5_insomnia_full.FUMA.txt.gz"
RSID = re.compile(r"^rs\d+$", re.IGNORECASE)
DNA = {"A", "C", "G", "T"}
MIN_P = 1e-300
N_CASES, N_CONTROLS = 51643, 446273


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def canonical_chrom(value: str) -> str | None:
    value = value.removeprefix("chr")
    return value if value.isdigit() and 1 <= int(value) <= 22 else None


def read_chain(path: Path) -> dict[str, tuple[list[int], list[tuple[int, int, str, int, int, str]]]]:
    """Le cadeia UCSC hg38->hg19 em blocos indexados pela coordenada GRCh38."""
    blocks: dict[str, list[tuple[int, int, str, int, int, str]]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="ascii") as handle:
        header = None
        tpos = qpos = 0
        for line in handle:
            line = line.strip()
            if not line:
                header = None
                continue
            if line.startswith("chain "):
                x = line.split()
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                if len(x) != 13:
                    raise ValueError("Cabecalho chain UCSC inesperado")
                chrom = canonical_chrom(x[2])
                qname, qsize, qstrand = x[7], int(x[8]), x[9]
                tpos, qpos = int(x[5]), int(x[10])
                header = (chrom, qname, qsize, qstrand)
                continue
            if header is None:
                continue
            x = [int(v) for v in line.split()]
            size = x[0]
            chrom, qname, qsize, qstrand = header
            dest = canonical_chrom(qname)
            if chrom is not None and dest is not None:
                blocks[chrom].append((tpos, tpos + size, dest, qpos, qsize, qstrand))
            tpos += size
            qpos += size
            if len(x) == 3:
                tpos += x[1]
                qpos += x[2]
    indexed = {}
    for chrom, values in blocks.items():
        values.sort(key=lambda x: x[0])
        indexed[chrom] = ([x[0] for x in values], values)
    return indexed


def lift_one(chain, chrom: str, pos_1based: int) -> tuple[str, int] | None:
    starts, blocks = chain[chrom]
    pos0 = pos_1based - 1
    i = bisect.bisect_right(starts, pos0) - 1
    if i < 0:
        return None
    start, end, dest, qpos, qsize, qstrand = blocks[i]
    if not start <= pos0 < end:
        return None
    offset = pos0 - start
    if qstrand == "+":
        return dest, qpos + offset + 1
    # qpos e dado na sequencia reversa; converte de volta a coordenada plus-strand.
    return dest, qsize - (qpos + offset)


def first_rsid(value: str) -> str | None:
    for candidate in re.split(r"[,;]", value):
        candidate = candidate.strip()
        if RSID.fullmatch(candidate):
            return candidate
    return None


def main() -> None:
    if not SOURCE.exists():
        raise FileNotFoundError(f"Arquivo FinnGen R13 nao encontrado: {SOURCE}")
    if not CHAIN.exists():
        raise FileNotFoundError(f"Cadeia UCSC GRCh38->hg19 nao encontrada: {CHAIN}")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    chain = read_chain(CHAIN)
    report = {
        "study_id": "finngen_r13_f5_insomnia",
        "label": "FinnGen R13 F5_INSOMNIA",
        "input_file": str(SOURCE), "input_bytes": SOURCE.stat().st_size,
        "input_build": "GRCh38", "output_build": "GRCh37/hg19",
        "liftover_chain": str(CHAIN), "liftover_chain_bytes": CHAIN.stat().st_size,
        "sample_size": N_CASES + N_CONTROLS, "cases": N_CASES, "controls": N_CONTROLS,
        "started_at": utc_now(),
        "qc_policy": {
            "chromosomes": "autossomos 1-22", "alleles": "SNVs bialelicos A/C/G/T",
            "rsid": "primeiro rsID valido por variante; obrigatorio", "maf_min": 0.001,
            "pvalue": "0<P<=1; todas as variantes validas sao mantidas, requisito para MAGMA gene-based; P abaixo de 1e-300 e limitada a 1e-300",
            "build": "liftover UCSC GRCh38 para GRCh37/hg19 antes da exportacao",
            "info": "nao disponivel na sumstat publica; o FinnGen exclui MAC <10 na analise original",
        },
        "qc": {key: 0 for key in (
            "raw_records", "kept", "non_autosomal", "non_biallelic_or_non_snv", "missing_or_invalid_rsid",
            "invalid_statistics", "maf_below_0_001", "liftover_unmapped", "liftover_non_autosomal", "duplicates_adjacent",
            "p_capped_at_1e_300")},
    }
    partial = OUTDIR / ".finngen_r13_f5_insomnia.partial.txt.gz"
    last_key = None
    try:
        with gzip.open(SOURCE, "rt", encoding="utf-8") as source, gzip.open(partial, "wt", encoding="utf-8", newline="") as target:
            header = source.readline().rstrip("\n").lstrip("#").split("\t")
            required = ["chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"]
            if any(x not in header for x in required):
                raise ValueError(f"Colunas FinnGen ausentes: esperado {required}; recebido {header}")
            ix = {name: header.index(name) for name in required}
            target.write("SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\tN\n")
            for line in source:
                report["qc"]["raw_records"] += 1
                x = line.rstrip("\n").split("\t")
                try:
                    chrom = canonical_chrom(x[ix["chrom"]]); pos = int(x[ix["pos"]])
                except (IndexError, ValueError):
                    report["qc"]["non_autosomal"] += 1; continue
                if chrom is None or pos <= 0:
                    report["qc"]["non_autosomal"] += 1; continue
                ref, alt = x[ix["ref"]].upper(), x[ix["alt"]].upper()
                if ref not in DNA or alt not in DNA or ref == alt:
                    report["qc"]["non_biallelic_or_non_snv"] += 1; continue
                rsid = first_rsid(x[ix["rsids"]])
                if rsid is None:
                    report["qc"]["missing_or_invalid_rsid"] += 1; continue
                try:
                    p, beta, se, af = (float(x[ix[k]]) for k in ("pval", "beta", "sebeta", "af_alt"))
                except (ValueError, IndexError):
                    report["qc"]["invalid_statistics"] += 1; continue
                if not all(math.isfinite(v) for v in (p, beta, se, af)) or se <= 0 or not 0 <= af <= 1 or p < 0 or p > 1:
                    report["qc"]["invalid_statistics"] += 1; continue
                if min(af, 1 - af) < 0.001:
                    report["qc"]["maf_below_0_001"] += 1; continue
                if p < MIN_P:
                    p = MIN_P; report["qc"]["p_capped_at_1e_300"] += 1
                lifted = lift_one(chain, chrom, pos)
                if lifted is None:
                    report["qc"]["liftover_unmapped"] += 1; continue
                out_chrom, out_pos = lifted
                if out_chrom not in chain:
                    report["qc"]["liftover_non_autosomal"] += 1; continue
                key = f"{out_chrom}:{out_pos}:{alt}:{ref}"
                if key == last_key:
                    report["qc"]["duplicates_adjacent"] += 1; continue
                last_key = key
                target.write(f"{rsid}\t{out_chrom}\t{out_pos}\t{alt}\t{ref}\t{beta:.17g}\t{se:.17g}\t{p:.17g}\t{N_CASES + N_CONTROLS}\n")
                report["qc"]["kept"] += 1
        partial.replace(OUT)
    finally:
        partial.unlink(missing_ok=True)
    report["finished_at"] = utc_now(); report["fuma_file"] = str(OUT); report["fuma_file_bytes"] = OUT.stat().st_size
    report["output_columns"] = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"]
    (OUTDIR / "finngen_r13_f5_insomnia_full.qc.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    metadata = {"study_id": report["study_id"], "label": report["label"],
                "role": "Validacao externa principal: coorte FinnGen independente do UK Biobank; nao redefine sozinha o core16.",
                "phenotype": "Insomnia clinica por registros de saude: ICD-10 F51.0 ou G47.0; controles excluem F5_SLEEP.",
                "population": "Finnish/European", "build": "GRCh38 convertido para GRCh37/hg19 por UCSC liftOver",
                "cases": N_CASES, "controls": N_CONTROLS, "N": N_CASES + N_CONTROLS,
                "effect_allele": "A1=ALT (beta FinnGen relativo ao alelo alternativo)", "fuma_file": OUT.name,
                "columns": report["output_columns"], "status": "VALIDATE_BEFORE_UPLOAD"}
    (OUTDIR / "finngen_r13_f5_insomnia_full.submission_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"output": str(OUT), "qc": report["qc"]}, indent=2))


if __name__ == "__main__":
    main()
