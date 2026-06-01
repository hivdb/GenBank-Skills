#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import tempfile
from collections import defaultdict
from io import StringIO
from multiprocessing import Pool
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


DEFAULT_MATCHED_FASTA_REPORT = Path("outputs/NS5B_matched_fasta_files.txt")
FASTA_EXTENSIONS = {".fa", ".faa", ".fasta", ".fna", ".fas", ".ffn", ".frn", ".seq"}
ACCESSION_RE = re.compile(r"^ACCESSION\s+(\S+)", re.MULTILINE)
COMMENT_RE = re.compile(r"^COMMENT\s+(.+?)(?=^FEATURES\s)", re.MULTILINE | re.DOTALL)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read all FASTA files listed in an NS5B matched-FASTA report, collect GenBank "
            "accessions, extract matching GenBank records from a local archive, and write "
            "source-feature qualifiers plus raw structured-comment text to a CSV."
        )
    )
    parser.add_argument("--matched-fasta-report", default=str(DEFAULT_MATCHED_FASTA_REPORT))
    parser.add_argument("--genbank-dir", required=True, help="Directory containing local GenBank flatfiles (*.seq)")
    parser.add_argument("--output-csv", default="", help="Output CSV path")
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 1) - 1),
        help="Parallel worker count for processing extracted GenBank records",
    )
    parser.add_argument("--keep-extracted-dir", action="store_true")
    return parser.parse_args()


def script_temp_dir() -> Path:
    path = Path("temp") / "hcv-ns5b-build-workflow" / Path(__file__).stem
    path.mkdir(parents=True, exist_ok=True)
    return path


def flatten_value(value: Any) -> str:
    if isinstance(value, list):
        return " | ".join(str(item) for item in value)
    return str(value)


def source_feature(record) -> SeqFeature | None:
    for feature in record.features:
        if feature.type == "source":
            return feature
    return None


def refid_from_fasta_path(path: Path) -> str:
    match = re.match(r"(\d+)", path.name)
    return match.group(1) if match else ""


def read_matched_fasta_paths(report_path: Path) -> list[Path]:
    return [Path(line.strip()).expanduser() for line in report_path.read_text(encoding="utf-8").splitlines() if line.strip()]


def collect_accessions_from_fastas(fasta_paths: list[Path]) -> tuple[set[str], dict[str, set[str]], dict[str, set[str]]]:
    accessions: set[str] = set()
    fasta_files_by_accession: dict[str, set[str]] = defaultdict(set)
    refids_by_accession: dict[str, set[str]] = defaultdict(set)
    for fasta_path in fasta_paths:
        if not fasta_path.is_file():
            raise RuntimeError(f"FASTA file was not found: {fasta_path}")
        if fasta_path.suffix.lower() not in FASTA_EXTENSIONS:
            continue
        refid = refid_from_fasta_path(fasta_path)
        for record in SeqIO.parse(fasta_path, "fasta"):
            accession = str(record.id).strip()
            if not accession:
                continue
            accessions.add(accession)
            fasta_files_by_accession[accession].add(str(fasta_path.resolve()))
            if refid:
                refids_by_accession[accession].add(refid)
    return accessions, fasta_files_by_accession, refids_by_accession


def iter_raw_genbank_records(path: Path):
    chunks: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            chunks.append(line)
            if line.rstrip("\n") == "//":
                yield "".join(chunks)
                chunks = []
    if chunks:
        tail = "".join(chunks).strip()
        if tail:
            yield "".join(chunks)


def accession_from_raw_record(record_text: str) -> str:
    match = ACCESSION_RE.search(record_text)
    return match.group(1).strip() if match else ""


def extract_matching_genbank_records(genbank_dir: Path, target_accessions: set[str], extracted_dir: Path) -> dict[str, str]:
    gb_file_by_accession: dict[str, str] = {}
    for seq_path in sorted(genbank_dir.glob("*.seq")):
        for record_text in iter_raw_genbank_records(seq_path):
            accession = accession_from_raw_record(record_text)
            if accession not in target_accessions:
                continue
            out_path = extracted_dir / f"{accession}.gb"
            out_path.write_text(record_text if record_text.endswith("\n") else f"{record_text}\n", encoding="utf-8")
            gb_file_by_accession[accession] = str(seq_path.resolve())
    return gb_file_by_accession


def extract_raw_structured_comment(record_text: str) -> str:
    comment_match = COMMENT_RE.search(record_text)
    if not comment_match:
        return ""
    lines = comment_match.group(1).splitlines()
    normalized = [line[12:] if len(line) >= 12 else line for line in lines]
    blocks: list[str] = []
    current: list[str] = []
    in_block = False
    for line in normalized:
        text = line.rstrip()
        if text.startswith("##") and text.endswith("-START##"):
            if current:
                blocks.append("\n".join(current).strip())
                current = []
            in_block = True
            current.append(text)
            continue
        if in_block:
            current.append(text)
            if text.startswith("##") and text.endswith("-END##"):
                blocks.append("\n".join(current).strip())
                current = []
                in_block = False
    if current:
        blocks.append("\n".join(current).strip())
    return "\n\n".join(block for block in blocks if block)


def process_extracted_record(task: tuple[str, str, str, str]) -> dict[str, str]:
    accession, gb_path_text, _fasta_files_text, refids_text = task
    gb_path = Path(gb_path_text)
    record_text = gb_path.read_text(encoding="utf-8")
    record = SeqIO.read(StringIO(record_text), "genbank")
    source = source_feature(record)
    qualifiers = source.qualifiers if source else {}
    row: dict[str, str] = {
        "RefID": refids_text,
        "Accession": accession,
        "definition": str(record.description or ""),
        "source_feature_present": "yes" if source else "no",
        "StructuredComment": extract_raw_structured_comment(record_text),
    }
    for key, value in sorted(qualifiers.items()):
        row[f"source_{key}"] = flatten_value(value)
    return row


def build_rows(
    extracted_dir: Path,
    target_accessions: set[str],
    fasta_files_by_accession: dict[str, set[str]],
    refids_by_accession: dict[str, set[str]],
    gb_archive_file_by_accession: dict[str, str],
    workers: int,
) -> tuple[list[dict[str, str]], list[str]]:
    rows: list[dict[str, str]] = []
    found_accessions = set(gb_archive_file_by_accession)
    tasks: list[tuple[str, str, str, str]] = []
    for accession in sorted(found_accessions):
        tasks.append(
            (
                accession,
                str((extracted_dir / f"{accession}.gb").resolve()),
                " | ".join(sorted(fasta_files_by_accession.get(accession, set()))),
                " | ".join(sorted(refids_by_accession.get(accession, set()))),
            )
        )
    if tasks:
        if workers <= 1:
            rows = [process_extracted_record(task) for task in tasks]
        else:
            with Pool(processes=workers) as pool:
                rows = list(pool.imap_unordered(process_extracted_record, tasks))
    missing_accessions = sorted(target_accessions - found_accessions)
    for accession in missing_accessions:
        rows.append(
            {
                "RefID": " | ".join(sorted(refids_by_accession.get(accession, set()))),
                "Accession": accession,
                "definition": "",
                "source_feature_present": "no",
                "StructuredComment": "",
                "status": "genbank_record_not_found",
            }
        )
    return rows, missing_accessions


def ordered_fieldnames(rows: list[dict[str, str]]) -> list[str]:
    preferred = ["RefID", "Accession", "definition", "source_feature_present", "StructuredComment", "status"]
    discovered = sorted({key for row in rows for key in row.keys() if key not in preferred})
    return preferred + discovered


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    matched_fasta_report = Path(args.matched_fasta_report).expanduser()
    genbank_dir = Path(args.genbank_dir).expanduser()
    output_csv = Path(args.output_csv).expanduser() if args.output_csv else script_temp_dir() / "NS5B_SourceFeatures.csv"
    workers = max(1, args.workers)
    if not matched_fasta_report.is_file():
        raise RuntimeError(f"Matched FASTA report was not found: {matched_fasta_report}")
    if not genbank_dir.is_dir():
        raise RuntimeError(f"GenBank directory was not found: {genbank_dir}")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    if output_csv.exists():
        print({"output_csv": str(output_csv.resolve()), "status": "skipped_existing_output"})
        return 0
    fasta_paths = read_matched_fasta_paths(matched_fasta_report)
    target_accessions, fasta_files_by_accession, refids_by_accession = collect_accessions_from_fastas(fasta_paths)
    extracted_dir = Path(tempfile.mkdtemp(prefix="ns5b_sourcefeatures_", dir=script_temp_dir()))
    try:
        gb_archive_file_by_accession = extract_matching_genbank_records(genbank_dir, target_accessions, extracted_dir)
        rows, missing_accessions = build_rows(
            extracted_dir, target_accessions, fasta_files_by_accession, refids_by_accession, gb_archive_file_by_accession, workers
        )
        rows.sort(key=lambda row: (row["RefID"], row["Accession"]))
        fieldnames = ordered_fieldnames(rows)
        write_csv(output_csv, rows, fieldnames)
        print(
            {
                "output_csv": str(output_csv.resolve()),
                "matched_fasta_file_count": len(fasta_paths),
                "unique_accession_count": len(target_accessions),
                "extracted_genbank_record_count": len(gb_archive_file_by_accession),
                "row_count": len(rows),
                "missing_accession_count": len(missing_accessions),
                "workers": workers,
                "extracted_dir": str(extracted_dir),
            }
        )
    finally:
        if not args.keep_extracted_dir:
            shutil.rmtree(extracted_dir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
