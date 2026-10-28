#!/usr/bin/env python3
"""Standardize VEP/snpEff annotated VCF into downstream WGS tables."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path


MISSING = {"", ".", "NA", "N/A", "na", "null", "NULL"}
IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, "NA": 0}
LOF_TERMS = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
}
MISSENSE_TERMS = {"missense_variant", "protein_altering_variant"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotated-vcf", required=True)
    parser.add_argument("--genome-version", default="GRCh38")
    parser.add_argument("--max-variants", type=int, default=0)
    parser.add_argument("--rare-af-threshold", type=float, default=0.01)
    parser.add_argument("--outdir", default=".")
    return parser.parse_args()


def clean(value: object) -> str:
    if value is None:
        return "NA"
    text = str(value).strip()
    return "NA" if text in MISSING else text


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.endswith(".gz") else open(path, "r", encoding="utf-8", errors="replace")


def parse_info(info: str) -> dict[str, str]:
    out = {}
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value
        else:
            out[item] = "1"
    return out


def parse_csq_header(line: str) -> list[str]:
    match = re.search(r"Format: ([^\">]+)", line)
    if not match:
        return []
    return [field.strip() for field in match.group(1).split("|")]


def parse_ann_header(line: str) -> list[str]:
    match = re.search(r"'([^']+)'", line)
    if not match:
        return []
    return [field.strip().replace(" ", "_") for field in match.group(1).split("|")]


def split_annotation(value: str, fields: list[str]) -> list[dict[str, str]]:
    annotations = []
    for item in value.split(","):
        parts = item.split("|")
        row = {field: clean(parts[idx] if idx < len(parts) else "NA") for idx, field in enumerate(fields)}
        annotations.append(row)
    return annotations


def first_float(values: list[str]) -> str:
    """Return the first numeric value, including Number=A INFO fields."""
    for raw_value in values:
        for value in str(raw_value).split(","):
            value = value.strip()
            if value in MISSING:
                continue
            try:
                return f"{float(value):.6g}"
            except ValueError:
                continue
    return "NA"


def max_float(values: list[str]) -> str:
    parsed = []
    for value in values:
        if value in MISSING:
            continue
        try:
            parsed.append(float(value))
        except ValueError:
            continue
    return f"{max(parsed):.6g}" if parsed else "NA"


def first_present(mapping: dict[str, str], keys: list[str]) -> str:
    """Return the first non-missing value for a list of possible VEP names."""
    for key in keys:
        value = clean(mapping.get(key, "NA"))
        if value != "NA":
            return value
    return "NA"


def best_csq(csq_rows: list[dict[str, str]]) -> dict[str, str]:
    def score(row: dict[str, str]) -> tuple[int, int, int, int]:
        impact = IMPACT_RANK.get(row.get("IMPACT", "NA"), 0)
        mane = 1 if row.get("MANE_SELECT", "NA") != "NA" or row.get("MANE", "NA") != "NA" else 0
        canonical = 1 if row.get("CANONICAL", "NA") == "YES" else 0
        gene = 1 if row.get("SYMBOL", "NA") != "NA" else 0
        return impact, mane, canonical, gene

    return sorted(csq_rows, key=score, reverse=True)[0] if csq_rows else {}


def best_ann(ann_rows: list[dict[str, str]]) -> dict[str, str]:
    def score(row: dict[str, str]) -> tuple[int, int]:
        impact = IMPACT_RANK.get(row.get("Annotation_Impact", "NA"), 0)
        gene = 1 if row.get("Gene_Name", "NA") != "NA" else 0
        return impact, gene

    return sorted(ann_rows, key=score, reverse=True)[0] if ann_rows else {}


def normalize_record(chrom: str, pos: str, var_id: str, ref: str, alt: str, info: str, csq_fields: list[str], ann_fields: list[str], rare_af_threshold: float) -> dict[str, str]:
    info_map = parse_info(info)
    csq_rows = split_annotation(info_map.get("CSQ", ""), csq_fields) if info_map.get("CSQ") and csq_fields else []
    ann_rows = split_annotation(info_map.get("ANN", ""), ann_fields) if info_map.get("ANN") and ann_fields else []
    csq = best_csq(csq_rows)
    ann = best_ann(ann_rows)

    consequence = clean(csq.get("Consequence") or ann.get("Annotation"))
    impact = clean(csq.get("IMPACT") or ann.get("Annotation_Impact"))
    gene = clean(csq.get("SYMBOL") or ann.get("Gene_Name"))
    transcript = clean(csq.get("Feature") or ann.get("Feature_ID"))
    variant_id = var_id if clean(var_id) != "NA" else f"{chrom}:{pos}:{ref}:{alt}"

    af_values = [
        info_map.get("gnomADg_AF", "NA"),
        info_map.get("gnomAD_AF", "NA"),
        info_map.get("gnomad_AF", "NA"),
        info_map.get("AF", "NA"),
        csq.get("AF", "NA"),
        csq.get("EAS_AF", "NA"),
        csq.get("gnomADg_AF", "NA"),
        csq.get("gnomADg_AF_eas", "NA"),
    ]
    public_af = first_float(af_values)
    # These fields preserve the evidence required for formal gene masks. They
    # remain NA when the VEP run lacks the respective plugin/database.
    gnomad_eas_af = first_float(
        [
            info_map.get("gnomADg_AF_eas", "NA"),
            info_map.get("gnomADg_AF_EAS", "NA"),
            info_map.get("gnomAD_EAS_AF", "NA"),
            info_map.get("gnomad_EAS_AF", "NA"),
            csq.get("gnomADg_AF_eas", "NA"),
            csq.get("EAS_AF", "NA"),
        ]
    )
    loftee = first_present(csq, ["LoF", "LOFTEE", "LoFTEE"])
    revel = first_float([csq.get("REVEL", "NA"), csq.get("REVEL_score", "NA"), info_map.get("REVEL", "NA")])
    alphamissense = first_float(
        [
            csq.get("AlphaMissense", "NA"),
            csq.get("AlphaMissense_score", "NA"),
            csq.get("AlphaMissense_pred", "NA"),
            csq.get("am_pathogenicity", "NA"),
            info_map.get("AlphaMissense", "NA"),
        ]
    )
    alphamissense_class = first_present(csq, ["am_class", "AlphaMissense_class", "AlphaMissense_pred"])
    spliceai = max_float(
        [
            csq.get("SpliceAI_pred_DS_AG", "NA"),
            csq.get("SpliceAI_pred_DS_AL", "NA"),
            csq.get("SpliceAI_pred_DS_DG", "NA"),
            csq.get("SpliceAI_pred_DS_DL", "NA"),
        ]
    )
    cadd = first_float(
        [
            csq.get("CADD_phred", "NA"),
            csq.get("CADD_PHRED", "NA"),
            info_map.get("CADD", "NA"),
            info_map.get("CADD_phred", "NA"),
        ]
    )
    clinvar_clnsig = first_present(csq, ["ClinVar_CLNSIG", "CLNSIG", "ClinVar_CLIN_SIG"])
    clinvar_clnrevstat = first_present(csq, ["ClinVar_CLNREVSTAT", "CLNREVSTAT"])
    clinvar_clndn = first_present(csq, ["ClinVar_CLNDN", "CLNDN"])
    clinvar_clnsigconf = first_present(csq, ["ClinVar_CLNSIGCONF", "CLNSIGCONF"])
    gnomad_popmax_af = first_float(
        [
            csq.get("gnomADg_popmax_AF", "NA"),
            csq.get("gnomADg_AF_popmax", "NA"),
            info_map.get("gnomADg_popmax_AF", "NA"),
        ]
    )
    maxentscan = first_float([csq.get("MaxEntScan_diff", "NA"), csq.get("MaxEntScan", "NA"), info_map.get("MaxEntScan_diff", "NA")])
    spliceai_raw = first_present(
        csq,
        [
            "SpliceAI_pred",
            "SpliceAI",
            "SpliceAI_pred_DS_AG",
            "SpliceAI_pred_DS_AL",
            "SpliceAI_pred_DS_DG",
            "SpliceAI_pred_DS_DL",
        ],
    )
    sift = first_present(csq, ["SIFT", "SIFT_pred"])
    polyphen = first_present(csq, ["PolyPhen", "PolyPhen_pred"])
    terms = set(consequence.split("&")) if consequence != "NA" else set()
    is_rare = "NA"
    if public_af != "NA":
        is_rare = "1" if float(public_af) < rare_af_threshold else "0"

    return {
        "variant_id": variant_id,
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "gene": gene,
        "transcript": transcript,
        "consequence": consequence,
        "impact": impact,
        "cadd": cadd,
        "spliceai_delta": spliceai,
        "gnomad_af": public_af,
        "gnomad_eas_af": gnomad_eas_af,
        "loftee": loftee,
        "revel": revel,
        "alphamissense": alphamissense,
        "alphamissense_class": alphamissense_class,
        "clinvar_clnsig": clinvar_clnsig,
        "clinvar_clnrevstat": clinvar_clnrevstat,
        "clinvar_clndn": clinvar_clndn,
        "clinvar_clnsigconf": clinvar_clnsigconf,
        "gnomad_popmax_af": gnomad_popmax_af,
        "maxentscan": maxentscan,
        "spliceai_raw": spliceai_raw,
        "sift": sift,
        "polyphen": polyphen,
        "is_rare": is_rare,
        "is_lof": "1" if terms & LOF_TERMS else "0",
        "is_missense": "1" if terms & MISSENSE_TERMS else "0",
    }


def write_tsv(path: Path, header: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: clean(row.get(col, "NA")) for col in header})


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    results_dir = outdir / "results"
    reports_dir = outdir / "reports"
    logs_dir = outdir / "logs"
    for directory in (results_dir, reports_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    header = [
        "variant_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "gene",
        "transcript",
        "consequence",
        "impact",
        "cadd",
        "spliceai_delta",
        "gnomad_af",
        "gnomad_eas_af",
        "loftee",
        "revel",
        "alphamissense",
        "alphamissense_class",
        "clinvar_clnsig",
        "clinvar_clnrevstat",
        "clinvar_clndn",
        "clinvar_clnsigconf",
        "gnomad_popmax_af",
        "maxentscan",
        "spliceai_raw",
        "sift",
        "polyphen",
        "is_rare",
        "is_lof",
        "is_missense",
    ]
    warnings: list[dict[str, str]] = []
    stats = Counter()
    csq_fields: list[str] = []
    ann_fields: list[str] = []

    gene_map_header = [key for key in header if key not in {"chrom", "pos", "ref", "alt", "transcript"}]
    variant_path = results_dir / "variant_annotation.tsv"
    gene_path = results_dir / "gene_variant_map.tsv"
    row_count = 0
    with variant_path.open("w", encoding="utf-8", newline="") as variant_handle, gene_path.open(
        "w", encoding="utf-8", newline=""
    ) as gene_handle:
        variant_writer = csv.DictWriter(variant_handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        gene_writer = csv.DictWriter(gene_handle, fieldnames=gene_map_header, delimiter="\t", lineterminator="\n")
        variant_writer.writeheader()
        gene_writer.writeheader()
        with open_text(args.annotated_vcf) as handle:
            for line in handle:
                if line.startswith("##INFO=<ID=CSQ"):
                    csq_fields = parse_csq_header(line)
                    stats["has_csq_header"] = 1
                    continue
                if line.startswith("##INFO=<ID=ANN"):
                    ann_fields = parse_ann_header(line)
                    stats["has_ann_header"] = 1
                    continue
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    stats["records_skipped_malformed"] += 1
                    continue
                chrom, pos, var_id, ref, alt, _qual, _filt, info = parts[:8]
                row = normalize_record(chrom, pos, var_id, ref, alt, info, csq_fields, ann_fields, args.rare_af_threshold)
                variant_writer.writerow({col: clean(row.get(col, "NA")) for col in header})
                if row["gene"] != "NA":
                    gene_writer.writerow({col: clean(row.get(col, "NA")) for col in gene_map_header})
                row_count += 1
                stats["records_standardized"] += 1
                if row["gene"] == "NA":
                    stats["records_without_gene"] += 1
                if row["is_lof"] == "1":
                    stats["lof_variants"] += 1
                if row["is_missense"] == "1":
                    stats["missense_variants"] += 1
                if row["is_rare"] == "1":
                    stats["rare_variants"] += 1
                if row["loftee"] != "NA":
                    stats["variants_with_loftee"] += 1
                if row["spliceai_delta"] != "NA":
                    stats["variants_with_spliceai"] += 1
                if row["gnomad_eas_af"] != "NA":
                    stats["variants_with_gnomad_eas_af"] += 1
                if row["clinvar_clnsig"] != "NA":
                    stats["variants_with_clinvar"] += 1
                if row["alphamissense"] != "NA":
                    stats["variants_with_alphamissense"] += 1
                if row["maxentscan"] != "NA":
                    stats["variants_with_maxentscan"] += 1
                if row["cadd"] != "NA":
                    stats["variants_with_cadd"] += 1
                if row["revel"] != "NA":
                    stats["variants_with_revel"] += 1
                if args.max_variants and row_count >= args.max_variants:
                    break

    if not csq_fields and not ann_fields:
        warnings.append(
            {
                "level": "WARNING",
                "code": "NO_FUNCTIONAL_ANNOTATION_HEADER",
                "message": "No CSQ or ANN header found; output contains variant coordinates but annotation fields may be NA.",
            }
        )
    if row_count == 0:
        raise RuntimeError("No VCF records were standardized.")
    summary_rows = [{"metric": key, "value": value} for key, value in sorted(stats.items())]
    summary_rows.extend(
        [
            {"metric": "genome_version", "value": args.genome_version},
            {"metric": "max_variants", "value": args.max_variants},
            {"metric": "rare_af_threshold", "value": args.rare_af_threshold},
            {"metric": "result_rows", "value": row_count},
            {"metric": "gene_map_rows", "value": row_count - stats.get("records_without_gene", 0)},
        ]
    )
    write_tsv(reports_dir / "variant_annotation_summary.tsv", ["metric", "value"], summary_rows)
    write_tsv(logs_dir / "variant_annotation_warnings.tsv", ["level", "code", "message"], warnings)
    (reports_dir / "variant_annotation_summary.json").write_text(
        json.dumps({"status": "PASS", "stats": dict(stats), "result_rows": row_count, "warning_count": len(warnings)}, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    report = [
        "# WGS variant annotation standardization report",
        "",
        "Status: PASS",
        f"Annotated VCF: {args.annotated_vcf}",
        f"Genome version: {args.genome_version}",
        f"Variants standardized: {row_count}",
        f"Gene-map rows: {row_count - stats.get('records_without_gene', 0)}",
        "",
        "## Key Counts",
        "",
        "| metric | value |",
        "|---|---:|",
        f"| rare_variants | {stats.get('rare_variants', 0)} |",
        f"| lof_variants | {stats.get('lof_variants', 0)} |",
        f"| missense_variants | {stats.get('missense_variants', 0)} |",
        f"| records_without_gene | {stats.get('records_without_gene', 0)} |",
        "",
        "## Warnings",
        "",
    ]
    if warnings:
        for warning in warnings:
            report.append(f"- {warning['code']}: {warning['message']}")
    else:
        report.append("No warnings.")
    (reports_dir / "variant_annotation_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
