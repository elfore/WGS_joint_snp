#!/usr/bin/env python3
"""Discover gVCFs and normalize their VCF headers for the WGS SNP workflow."""

import argparse
import csv
import gzip
import re
import sys
from contextlib import contextmanager
from pathlib import Path


CONTIG_PREFIX = "##contig=<ID="
VCF_SUFFIXES = (".g.vcf.gz", ".g.vcf.bgz", ".g.vcf", ".vcf.gz", ".vcf.bgz", ".vcf")
CONTIG_RE = re.compile(r"^##contig=<ID=([^,>]+)(?:,length=([0-9]+))?.*>$")


class PipelineError(RuntimeError):
    pass


@contextmanager
def open_text(path):
    path = Path(path)
    if path.name.endswith(".gz") or path.name.endswith(".bgz"):
        handle = gzip.open(str(path), "rt", encoding="utf-8", newline="")
    else:
        handle = open(str(path), "r", encoding="utf-8", newline="")
    try:
        yield handle
    finally:
        handle.close()


def ensure_parent(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def write_tsv(path, fieldnames, rows):
    ensure_parent(path)
    with open(str(path), "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_lines(path, lines):
    ensure_parent(path)
    with open(str(path), "w", encoding="utf-8", newline="") as handle:
        for line in lines:
            handle.write(line.rstrip("\r\n") + "\n")


def truthy(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def sample_id_from_filename(name):
    for suffix in VCF_SUFFIXES:
        if name.endswith(suffix):
            sample_id = name[:-len(suffix)]
            if sample_id:
                return sample_id
    raise PipelineError("Cannot derive sample ID from VCF filename: {}".format(name))


def parse_contig_line(line):
    match = CONTIG_RE.match(line.rstrip("\r\n"))
    if not match:
        return None
    contig_id, length = match.groups()
    return contig_id, int(length) if length else None


def read_vcf_header(path):
    header = []
    chrom_line = None
    with open_text(path) as handle:
        for line in handle:
            clean = line.rstrip("\r\n")
            if not clean.startswith("#"):
                break
            header.append(clean)
            if clean.startswith("#CHROM"):
                chrom_line = clean
                break
    if chrom_line is None:
        raise PipelineError("VCF header has no #CHROM line: {}".format(path))
    return header, chrom_line


def header_contigs(header):
    contigs = {}
    order = []
    for line in header:
        parsed = parse_contig_line(line)
        if parsed is None:
            continue
        contig_id, length = parsed
        if contig_id not in contigs:
            order.append(contig_id)
        contigs[contig_id] = length
    return order, contigs


def header_samples(chrom_line):
    fields = chrom_line.split("\t")
    if len(fields) < 9:
        raise PipelineError("Malformed #CHROM line: {}".format(chrom_line))
    return fields[9:]


def manifest_rows(path):
    with open(str(path), "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "sample_id" not in reader.fieldnames or "gvcf" not in reader.fieldnames:
            raise PipelineError("Manifest must contain sample_id and gvcf columns: {}".format(path))
        rows = []
        seen = set()
        for row in reader:
            sample_id = (row.get("sample_id") or "").strip()
            gvcf = (row.get("gvcf") or "").strip()
            if not sample_id or not gvcf:
                raise PipelineError("Manifest contains an empty sample_id or gvcf: {}".format(path))
            if sample_id in seen:
                raise PipelineError("Duplicate sample_id in manifest: {}".format(sample_id))
            seen.add(sample_id)
            rows.append({"sample_id": sample_id, "gvcf": gvcf})
    if not rows:
        raise PipelineError("Manifest contains no gVCF records: {}".format(path))
    return rows


def discover_gvcfs(args):
    root = Path(args.gvcf_dir)
    if not root.exists():
        raise PipelineError("gVCF directory does not exist: {}".format(root))
    if not root.is_dir():
        raise PipelineError("gVCF path is not a directory: {}".format(root))

    files = []
    for path in root.iterdir():
        if path.name.startswith(".") or not path.is_file():
            continue
        if path.name.endswith(VCF_SUFFIXES):
            files.append(path)
    files.sort(key=lambda item: item.name)
    if not files:
        raise PipelineError("No VCF/gVCF files found in {}".format(root))

    rows = []
    seen = set()
    for path in files:
        sample_id = sample_id_from_filename(path.name)
        if sample_id in seen:
            raise PipelineError("Duplicate sample ID derived from filenames: {}".format(sample_id))
        seen.add(sample_id)
        rows.append({"sample_id": sample_id, "gvcf": str(path.absolute())})

    write_tsv(args.out_manifest, ["sample_id", "gvcf"], rows)
    print("Discovered {} gVCF files".format(len(rows)))


def read_fai(path):
    contigs = []
    lengths = {}
    with open(str(path), "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\r\n").split("\t")
            if not fields or not fields[0]:
                continue
            if len(fields) < 2:
                raise PipelineError("Malformed FAI line {} in {}".format(line_number, path))
            contig_id = fields[0]
            try:
                length = int(fields[1])
            except ValueError:
                raise PipelineError("Invalid FAI length at {}:{}".format(path, line_number))
            if contig_id not in lengths:
                contigs.append(contig_id)
            lengths[contig_id] = length
    if not contigs:
        raise PipelineError("Reference FAI is empty: {}".format(path))
    return contigs, lengths


def resolve_reference_fai(reference_fai, reference_fasta, auto_detect):
    candidates = []
    if reference_fai:
        candidates.append(Path(reference_fai))
    if reference_fasta:
        fasta = Path(reference_fasta)
        candidates.append(Path(str(fasta) + ".fai"))
    for candidate in candidates:
        if candidate.is_file():
            return candidate, "found"
    if auto_detect:
        searched = ", ".join(str(item) for item in candidates) or "<none>"
        raise PipelineError("Unable to auto-detect a reference FAI; checked {}".format(searched))
    return None, "not_found"


def collect_input_metadata(rows):
    metadata = []
    union_order = []
    union_seen = set()
    common = None
    for row in rows:
        path = Path(row["gvcf"])
        if not path.is_file():
            raise PipelineError("gVCF listed in manifest does not exist: {}".format(path))
        header, chrom_line = read_vcf_header(path)
        order, lengths = header_contigs(header)
        if not order:
            raise PipelineError("gVCF has no ##contig lines: {}".format(path))
        samples = header_samples(chrom_line)
        for contig_id in order:
            if contig_id not in union_seen:
                union_seen.add(contig_id)
                union_order.append(contig_id)
        contig_set = set(order)
        common = contig_set if common is None else common.intersection(contig_set)
        metadata.append({
            "sample_id": row["sample_id"],
            "gvcf": str(path),
            "header": header,
            "chrom_line": chrom_line,
            "contig_order": order,
            "contig_lengths": lengths,
            "samples": samples,
        })
    return metadata, union_order, common or set()


def choose_template_contigs(mode, reference_order, reference_lengths, input_order,
                            input_common, input_union):
    mode = mode.strip().lower()
    if mode in {"reference_intersection", "ref_intersection"}:
        if not reference_order:
            raise PipelineError("reference_intersection requires a reference FAI")
        return [item for item in reference_order if item in input_union]
    if mode in {"reference", "reference_all"}:
        if not reference_order:
            raise PipelineError("{} requires a reference FAI".format(mode))
        return list(reference_order)
    if mode in {"input_intersection", "intersection"}:
        return [item for item in (reference_order or input_order) if item in input_common]
    if mode in {"input_union", "union"}:
        if reference_order:
            return [item for item in reference_order if item in input_union]
        return list(input_order)
    raise PipelineError("Unsupported contig mode: {}".format(mode))


def build_template(args):
    rows = manifest_rows(args.manifest)
    metadata, input_order, input_common = collect_input_metadata(rows)
    input_union = set(input_order)

    fai_path, fai_status = resolve_reference_fai(
        args.reference_fai, args.reference_fasta, truthy(args.auto_detect_reference_fai)
    )
    reference_order = []
    reference_lengths = {}
    if fai_path is not None:
        reference_order, reference_lengths = read_fai(fai_path)

    template_order = choose_template_contigs(
        args.contig_mode, reference_order, reference_lengths,
        input_order, input_common, input_union
    )
    if not template_order:
        raise PipelineError("No contigs remain after applying contig mode {}".format(args.contig_mode))

    input_length_maps = [item["contig_lengths"] for item in metadata]
    template_lines = []
    template_rows = []
    for index, contig_id in enumerate(template_order, 1):
        length = reference_lengths.get(contig_id)
        if length is None:
            for length_map in input_length_maps:
                if length_map.get(contig_id) is not None:
                    length = length_map[contig_id]
                    break
        if length is None:
            raise PipelineError("No length available for template contig {}".format(contig_id))
        template_lines.append("##contig=<ID={},length={}>".format(contig_id, length))
        template_rows.append({
            "order": index,
            "contig": contig_id,
            "length": length,
            "source": "reference_fai" if contig_id in reference_lengths else "input_header",
        })

    template_set = set(template_order)
    audit_rows = []
    for item in metadata:
        input_set = set(item["contig_order"])
        missing = [contig for contig in template_order if contig not in input_set]
        extra = [contig for contig in item["contig_order"] if contig not in template_set]
        mismatches = []
        for contig_id in template_order:
            observed = item["contig_lengths"].get(contig_id)
            expected = reference_lengths.get(contig_id)
            if observed is not None and expected is not None and observed != expected:
                mismatches.append("{}:{}!={}".format(contig_id, observed, expected))
        audit_rows.append({
            "sample_id": item["sample_id"],
            "gvcf": item["gvcf"],
            "input_contig_count": len(item["contig_order"]),
            "input_contigs": ",".join(item["contig_order"]),
            "missing_template_contigs": ",".join(missing),
            "extra_input_contigs": ",".join(extra),
            "length_mismatches": ",".join(mismatches),
            "header_samples": ",".join(item["samples"]),
        })

    excluded_rows = []
    for contig_id in reference_order:
        if contig_id not in input_union:
            excluded_rows.append({
                "contig": contig_id,
                "length": reference_lengths[contig_id],
                "reason": "absent_from_input_headers",
            })
    for contig_id in input_order:
        if contig_id not in set(reference_order) and reference_order:
            observed_length = ""
            for item in metadata:
                if contig_id in item["contig_lengths"]:
                    observed_length = item["contig_lengths"][contig_id]
                    break
            excluded_rows.append({
                "contig": contig_id,
                "length": observed_length,
                "reason": "not_in_reference_fai",
            })

    resolution_rows = [{
        "reference_fai": str(fai_path) if fai_path else "",
        "reference_fasta": args.reference_fasta,
        "reference_fai_status": fai_status,
        "auto_detect_reference_fai": str(truthy(args.auto_detect_reference_fai)).lower(),
        "contig_mode": args.contig_mode,
        "input_file_count": len(metadata),
        "input_union_count": len(input_union),
        "input_intersection_count": len(input_common),
        "template_contig_count": len(template_order),
    }]

    write_lines(args.out_contigs, template_lines)
    write_tsv(args.out_template_tsv, ["order", "contig", "length", "source"], template_rows)
    write_tsv(
        args.out_audit_tsv,
        ["sample_id", "gvcf", "input_contig_count", "input_contigs",
         "missing_template_contigs", "extra_input_contigs",
         "length_mismatches", "header_samples"],
        audit_rows,
    )
    write_tsv(args.out_excluded_tsv, ["contig", "length", "reason"], excluded_rows)
    write_tsv(
        args.out_reference_tsv,
        ["reference_fai", "reference_fasta", "reference_fai_status",
         "auto_detect_reference_fai", "contig_mode", "input_file_count",
         "input_union_count", "input_intersection_count", "template_contig_count"],
        resolution_rows,
    )
    print("Built template with {} contigs from {} gVCFs".format(len(template_order), len(metadata)))


def read_template_contigs(path):
    lines = []
    contig_ids = []
    seen = set()
    with open(str(path), "r", encoding="utf-8") as handle:
        for line in handle:
            clean = line.rstrip("\r\n")
            if not clean:
                continue
            parsed = parse_contig_line(clean)
            if parsed is None:
                raise PipelineError("Template contains a non-contig line: {}".format(clean))
            contig_id, _ = parsed
            if contig_id in seen:
                raise PipelineError("Template contains duplicate contig: {}".format(contig_id))
            seen.add(contig_id)
            contig_ids.append(contig_id)
            lines.append(clean)
    if not lines:
        raise PipelineError("Contig template is empty: {}".format(path))
    return lines, contig_ids


def rewrite_header(args):
    template_lines, template_order = read_template_contigs(args.template_contigs)
    header, chrom_line = read_vcf_header(args.gvcf)
    original_order, original_lengths = header_contigs(header)
    original_samples = header_samples(chrom_line)
    if not original_samples:
        raise PipelineError("Input gVCF has no sample column: {}".format(args.gvcf))

    rewritten = []
    inserted = False
    for line in header:
        if parse_contig_line(line) is not None:
            if not inserted:
                rewritten.extend(template_lines)
                inserted = True
            continue
        if line.startswith("#CHROM"):
            if not inserted:
                rewritten.extend(template_lines)
                inserted = True
            fields = line.split("\t")
            fields[9:] = [args.sample_id] * len(fields[9:])
            rewritten.append("\t".join(fields))
        else:
            rewritten.append(line)
    if not inserted:
        raise PipelineError("Could not insert contig template into VCF header")

    write_lines(args.out_header, rewritten)
    template_set = set(template_order)
    original_set = set(original_order)
    added = [item for item in template_order if item not in original_set]
    removed = [item for item in original_order if item not in template_set]
    length_mismatches = []
    for contig_id in template_order:
        if contig_id in original_lengths:
            template_length = parse_contig_line(
                next(line for line in template_lines
                     if parse_contig_line(line)[0] == contig_id)
            )[1]
            if original_lengths[contig_id] is not None and template_length is not None:
                if original_lengths[contig_id] != template_length:
                    length_mismatches.append(
                        "{}:{}!={}".format(contig_id, original_lengths[contig_id], template_length)
                    )

    report_row = {
        "sample_id": args.sample_id,
        "gvcf": str(args.gvcf),
        "original_sample_names": ",".join(original_samples),
        "input_contig_count": len(original_order),
        "template_contig_count": len(template_order),
        "added_contigs": ",".join(added),
        "removed_contigs": ",".join(removed),
        "length_mismatches": ",".join(length_mismatches),
        "status": "ok",
    }
    write_tsv(
        args.report_tsv,
        ["sample_id", "gvcf", "original_sample_names", "input_contig_count",
         "template_contig_count", "added_contigs", "removed_contigs",
         "length_mismatches", "status"],
        [report_row],
    )
    print("Rewrote header for {}".format(args.sample_id))


def build_parser():
    parser = argparse.ArgumentParser(
        description="Utilities for discovering gVCFs and normalizing VCF headers."
    )
    subparsers = parser.add_subparsers(dest="command")

    discover = subparsers.add_parser("discover-gvcfs")
    discover.add_argument("--gvcf-dir", required=True)
    discover.add_argument("--out-manifest", required=True)
    discover.set_defaults(handler=discover_gvcfs)

    template = subparsers.add_parser("build-template")
    template.add_argument("--manifest", required=True)
    template.add_argument("--reference-fai", required=True)
    template.add_argument("--reference-fasta", required=True)
    template.add_argument("--auto-detect-reference-fai", required=True)
    template.add_argument("--contig-mode", required=True)
    template.add_argument("--out-contigs", required=True)
    template.add_argument("--out-template-tsv", required=True)
    template.add_argument("--out-audit-tsv", required=True)
    template.add_argument("--out-excluded-tsv", required=True)
    template.add_argument("--out-reference-tsv", required=True)
    template.set_defaults(handler=build_template)

    rewrite = subparsers.add_parser("rewrite-header")
    rewrite.add_argument("--sample-id", required=True)
    rewrite.add_argument("--gvcf", required=True)
    rewrite.add_argument("--template-contigs", required=True)
    rewrite.add_argument("--out-header", required=True)
    rewrite.add_argument("--report-tsv", required=True)
    rewrite.set_defaults(handler=rewrite_header)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "handler"):
        parser.print_help(sys.stderr)
        return 2
    try:
        args.handler(args)
    except (OSError, csv.Error, PipelineError) as exc:
        print("ERROR: {}".format(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
