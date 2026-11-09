"""Converte VCFs OpenGWAS em inputs FUMA validados por QC.

Os VCFs são lidos por streaming: não são carregados por inteiro na memória.
O alelo ALT do VCF é sempre exportado como A1 (alelo de efeito de ES).
Uso:
  python projeto_insonia_tcc/pipeline_fuma/process_opengwas_vcf_to_fuma.py --study finn_b_f5_insomnia
  python projeto_insonia_tcc/pipeline_fuma/process_opengwas_vcf_to_fuma.py --study all
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import sqlite3
import tempfile
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = ROOT / "projeto_insonia_tcc" / "fuma_inputs_opengwas"
RSID = re.compile(r"^rs\d+$", re.IGNORECASE)
DNA = {"A", "C", "G", "T"}
MIN_P = 1e-300  # evita P=0, que o FUMA não aceita; só afeta valores abaixo desse limite numérico.


@dataclass(frozen=True)
class Study:
    key: str
    vcf: Path
    label: str
    phenotype: str
    population: str
    build: str
    fixed_n: int | None
    fixed_cases: int | None
    fixed_controls: int | None
    effect_note: str
    role: str


STUDIES = {
    "finn_b_f5_insomnia": Study(
        key="finn_b_f5_insomnia",
        vcf=ROOT / "finn-b-F5_INSOMNIA.vcf" / "finn-b-F5_INSOMNIA.vcf",
        label="FinnGen F5_INSOMNIA",
        phenotype="Insônia clínica registrada (endpoint FinnGen F5)",
        population="Finnish/European",
        build="GRCh37",
        fixed_n=392423,
        fixed_cases=4214,
        fixed_controls=388209,
        effect_note=("ES é exportado como BETA relativo ao ALT, conforme o VCF OpenGWAS. "
                     "A escala do efeito deve ser interpretada com cautela: o VCF rotula o StudyType como Continuous, "
                     "embora o endpoint seja clínico; esta análise é validação externa baseada em P-values/loci."),
        role="Validação externa; não é usado para redefinir a lista consenso principal.",
    ),
    "ukb_b_3957": Study(
        key="ukb_b_3957",
        vcf=ROOT / "ukb-b-3957 (1).vcf" / "ukb-b-3957 (1).vcf",
        label="OpenGWAS ukb-b-3957",
        phenotype="Sleeplessness / insomnia (UK Biobank, PHESANT)",
        population="European",
        build="GRCh37",
        fixed_n=462341,
        fixed_cases=None,
        fixed_controls=None,
        effect_note="ES é exportado como BETA relativo ao ALT; fenótipo UKB ordinal padronizado, não interpretar como log(OR).",
        role="Análise de sensibilidade UK Biobank; não representa evidência independente adicional.",
    ),
    "ukb_a_13": Study(
        key="ukb_a_13",
        vcf=ROOT / "ukb-a-13.vcf" / "ukb-a-13.vcf",
        label="OpenGWAS ukb-a-13 (Neale Lab)",
        phenotype="Sleeplessness / insomnia (UK Biobank)",
        population="European",
        build="GRCh37",
        fixed_n=336965,
        fixed_cases=None,
        fixed_controls=None,
        effect_note="ES é exportado como BETA relativo ao ALT; fenótipo UKB ordinal padronizado, não interpretar como log(OR).",
        role="Análise de sensibilidade UK Biobank; não representa evidência independente adicional.",
    ),
}


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def parse_float(value: str | None) -> float | None:
    if value in (None, "", ".", "NA"):
        return None
    try:
        x = float(value)
    except ValueError:
        return None
    return x if math.isfinite(x) else None


def parse_int(value: str | None) -> int | None:
    x = parse_float(value)
    if x is None or x <= 0:
        return None
    return int(round(x))


def vcf_header(path: Path) -> tuple[list[str], str, str]:
    sample_meta = ""
    column_header = ""
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##SAMPLE="):
                sample_meta = line.rstrip()
            if line.startswith("#CHROM"):
                column_header = line.rstrip()
                break
    if not column_header:
        raise ValueError(f"VCF sem cabeçalho #CHROM: {path}")
    return column_header.split("\t"), sample_meta, column_header


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def process(study: Study) -> dict:
    if not study.vcf.exists():
        raise FileNotFoundError(study.vcf)
    columns, sample_meta, _ = vcf_header(study.vcf)
    if len(columns) < 10:
        raise ValueError("VCF precisa ter ao menos 10 colunas")

    out_dir = OUTPUT_ROOT / study.key
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{study.key}.FUMA.txt.gz"
    partial_path = out_dir / f".{study.key}.partial.txt.gz"
    report = {
        "study_id": study.key, "label": study.label, "input_vcf": str(study.vcf),
        "input_bytes": study.vcf.stat().st_size, "build": study.build, "vcf_sample_metadata": sample_meta,
        "qc": {"raw_records": 0, "kept": 0, "non_autosomal": 0, "not_pass": 0, "non_biallelic_or_non_snv": 0,
               "invalid_rsid": 0, "missing_or_invalid_format": 0, "invalid_statistics": 0,
               "maf_below_0_001": 0, "info_below_0_8": 0, "duplicates_removed": 0, "p_capped_at_1e_300": 0},
        "started_at": now(), "qc_policy": {"chromosomes": "1-22", "filter": "PASS", "alleles": "biallelic SNVs only",
                                               "maf_min": 0.001, "info_min": "0.8 only when SI is present in record",
                                               "rsid": "required and validated as rs[digits]", "deduplication": "duplicatas adjacentes removidas; validação posterior verifica duplicatas globais"},
    }
    last_key = None
    try:
        with study.vcf.open("rt", encoding="utf-8") as source, gzip.open(partial_path, "wt", encoding="utf-8", newline="") as target:
            target.write("SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\tN\n")
            for line in source:
                if line.startswith("#"):
                    continue
                report["qc"]["raw_records"] += 1
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    report["qc"]["missing_or_invalid_format"] += 1; continue
                chrom, pos, rsid, ref, alt, _, filt, info, fmt, sample = fields[:10]
                try:
                    chrom_i, pos_i = int(chrom), int(pos)
                except ValueError:
                    report["qc"]["non_autosomal"] += 1; continue
                if chrom_i not in range(1, 23) or pos_i <= 0:
                    report["qc"]["non_autosomal"] += 1; continue
                if filt != "PASS":
                    report["qc"]["not_pass"] += 1; continue
                ref, alt = ref.upper(), alt.upper()
                if "," in alt or ref not in DNA or alt not in DNA:
                    report["qc"]["non_biallelic_or_non_snv"] += 1; continue
                if not RSID.fullmatch(rsid):
                    report["qc"]["invalid_rsid"] += 1; continue
                fmt_keys, values = fmt.split(":"), sample.split(":")
                if len(fmt_keys) != len(values):
                    report["qc"]["missing_or_invalid_format"] += 1; continue
                dat = dict(zip(fmt_keys, values))
                beta, se, lp, af = (parse_float(dat.get(x)) for x in ("ES", "SE", "LP", "AF"))
                n = parse_int(dat.get("SS")) or study.fixed_n
                if beta is None or se is None or se <= 0 or lp is None or lp < 0 or af is None or n is None:
                    report["qc"]["invalid_statistics"] += 1; continue
                maf = min(af, 1 - af)
                if maf < 0.001:
                    report["qc"]["maf_below_0_001"] += 1; continue
                si = parse_float(dat.get("SI"))
                if si is not None and si < 0.8:
                    report["qc"]["info_below_0_8"] += 1; continue
                p = 10.0 ** (-lp)
                if p < MIN_P:
                    p = MIN_P; report["qc"]["p_capped_at_1e_300"] += 1
                if not (0 < p <= 1):
                    report["qc"]["invalid_statistics"] += 1; continue
                key = f"{chrom_i}:{pos_i}:{alt}:{ref}"
                if key == last_key:
                    report["qc"]["duplicates_removed"] += 1; continue
                last_key = key
                target.write(f"{rsid}\t{chrom_i}\t{pos_i}\t{alt}\t{ref}\t{beta:.17g}\t{se:.17g}\t{p:.17g}\t{n}\n")
                report["qc"]["kept"] += 1
        partial_path.replace(out_path)
    finally:
        partial_path.unlink(missing_ok=True)
    report["finished_at"] = now(); report["fuma_candidate_file"] = str(out_path); report["fuma_candidate_bytes"] = out_path.stat().st_size
    report["output_columns"] = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"]
    write_json(out_dir / f"{study.key}.qc.json", report)
    metadata = {"study_id": study.key, "label": study.label, "role": study.role, "phenotype": study.phenotype,
                "population": study.population, "build": study.build, "fuma_file": out_path.name,
                "effect_allele": "A1 = ALT", "non_effect_allele": "A2 = REF", "effect_note": study.effect_note,
                "fixed_n_when_SS_absent": study.fixed_n, "cases": study.fixed_cases, "controls": study.fixed_controls,
                "columns": report["output_columns"], "status": "VALIDATE_BEFORE_UPLOAD"}
    write_json(out_dir / f"{study.key}.submission_metadata.json", metadata)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", choices=("all", *STUDIES), default="all")
    args = parser.parse_args()
    selected = STUDIES.values() if args.study == "all" else (STUDIES[args.study],)
    for study in selected:
        result = process(study)
        print(study.key, "DONE", result["qc"])


if __name__ == "__main__":
    main()
