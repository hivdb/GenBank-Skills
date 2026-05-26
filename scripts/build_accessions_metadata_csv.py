#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import tempfile
from io import StringIO
from multiprocessing import Pool
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

try:
    import tomllib
except ModuleNotFoundError as exc:  # pragma: no cover
    raise SystemExit(f"tomllib is required: {exc}") from exc


FASTA_EXTENSIONS = {".fa", ".faa", ".fasta", ".fna", ".fas", ".ffn", ".frn", ".seq"}
ACCESSION_RE = re.compile(r"^ACCESSION\s+(\S+)", re.MULTILINE)
COMMENT_RE = re.compile(r"^COMMENT\s+(.+?)(?=^FEATURES\s)", re.MULTILINE | re.DOTALL)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read all FASTA files in a directory, collect accession records across all FASTA files, "
            "extract matching GenBank records from a local archive, and write source-feature "
            "qualifiers plus raw structured-comment text to a CSV."
        )
    )
    parser.add_argument(
        "--fasta-dir",
        default="",
        help="Directory containing FASTA files with one or more accession records per file",
    )
    parser.add_argument(
        "--genbank-dir",
        default="",
        help="Directory containing local GenBank flatfiles (*.seq)",
    )
    parser.add_argument(
        "--pipeline-name",
        default="ns3",
        help="Pipeline section in pipeline.local.toml to use for default fasta/genbank paths",
    )
    parser.add_argument(
        "--output-csv",
        default="",
        help="Output CSV path",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 1) - 1),
        help="Parallel worker count for processing extracted GenBank records",
    )
    parser.add_argument(
        "--keep-extracted-dir",
        action="store_true",
        help="Keep the temporary directory containing extracted matching GenBank records",
    )
    return parser.parse_args()


def script_temp_dir() -> Path:
    path = Path("temp") / Path(__file__).stem
    path.mkdir(parents=True, exist_ok=True)
    return path


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def resolve_config_path(value: str) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return repo_root() / path


def load_pipeline_defaults(pipeline_name: str) -> dict[str, str]:
    config_path = repo_root() / "pipeline.local.toml"
    if not config_path.is_file():
        return {}

    with config_path.open("rb") as handle:
        data = tomllib.load(handle)

    merged: dict[str, object] = {}
    for section_name in ("common", pipeline_name):
        section = data.get(section_name, {})
        if isinstance(section, dict):
            merged.update(section)

    resolved: dict[str, str] = {}
    for key in ("fasta_pool", "genbank_dir"):
        value = merged.get(key)
        if value:
            resolved[key] = str(resolve_config_path(str(value)))
    return resolved


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
    return path.name.split("_", 1)[0] if "_" in path.name else ""


def iter_fasta_paths(fasta_dir: Path) -> list[Path]:
    return sorted(
        path for path in fasta_dir.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    )


def collect_accessions_from_fastas(fasta_paths: list[Path]) -> tuple[list[str], dict[str, str], dict[str, str], int]:
    accession_by_file: dict[str, str] = {}
    refid_by_accession: dict[str, str] = {}
    total_accession_records = 0

    for fasta_path in fasta_paths:
        refid = refid_from_fasta_path(fasta_path)
        for record in SeqIO.parse(fasta_path, "fasta"):
            accession = str(record.id).strip()
            if not accession:
                continue
            total_accession_records += 1
            if accession in accession_by_file:
                raise RuntimeError(f"Duplicate accession across FASTA files or records: {accession}")
            accession_by_file[accession] = str(fasta_path.resolve())
            refid_by_accession[accession] = refid

    accessions = sorted(accession_by_file)
    if not accessions:
        raise RuntimeError("No accession records were found across the FASTA files")
    if len(accessions) != total_accession_records:
        raise RuntimeError(
            f"Unique accession count ({len(accessions)}) does not match total FASTA accession record count ({total_accession_records})"
        )
    return accessions, accession_by_file, refid_by_accession, total_accession_records


def iter_raw_genbank_records(path: Path) -> Any:
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


def extract_matching_genbank_records(
    genbank_dir: Path,
    target_accessions: set[str],
    extracted_dir: Path,
) -> dict[str, str]:
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


def process_extracted_record(task: tuple[str, str, str]) -> dict[str, str]:
    accession, gb_path_text, refid = task
    gb_path = Path(gb_path_text)
    record_text = gb_path.read_text(encoding="utf-8")
    record = SeqIO.read(StringIO(record_text), "genbank")
    source = source_feature(record)
    qualifiers = source.qualifiers if source else {}

    row: dict[str, str] = {
        "RefID": refid,
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
    target_accessions: list[str],
    refid_by_accession: dict[str, str],
    gb_archive_file_by_accession: dict[str, str],
    workers: int,
) -> tuple[list[dict[str, str]], list[str]]:
    rows: list[dict[str, str]] = []
    found_accessions = set(gb_archive_file_by_accession)
    tasks: list[tuple[str, str, str]] = []

    for accession in sorted(found_accessions):
        tasks.append(
            (
                accession,
                str((extracted_dir / f"{accession}.gb").resolve()),
                refid_by_accession.get(accession, ""),
            )
        )

    if tasks:
        if workers <= 1:
            rows = [process_extracted_record(task) for task in tasks]
        else:
            with Pool(processes=workers) as pool:
                rows = list(pool.imap_unordered(process_extracted_record, tasks))

    missing_accessions = sorted(set(target_accessions) - found_accessions)
    for accession in missing_accessions:
        rows.append(
            {
                "RefID": refid_by_accession.get(accession, ""),
                "Accession": accession,
                "definition": "",
                "source_feature_present": "no",
                "StructuredComment": "",
                "status": "genbank_record_not_found",
            }
        )

    return rows, missing_accessions


def ordered_fieldnames(rows: list[dict[str, str]]) -> list[str]:
    preferred = [
        "RefID",
        "Accession",
        "definition",
        "source_feature_present",
        "StructuredComment",
        "status",
    ]
    discovered = sorted({key for row in rows for key in row.keys() if key not in preferred})
    return preferred + discovered


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    defaults = load_pipeline_defaults(args.pipeline_name)
    fasta_dir_value = args.fasta_dir or defaults.get("fasta_pool", "")
    genbank_dir_value = args.genbank_dir or defaults.get("genbank_dir", "")
    if not fasta_dir_value:
        raise RuntimeError(
            f"FASTA directory was not provided and no fasta_pool was found for pipeline '{args.pipeline_name}'"
        )
    if not genbank_dir_value:
        raise RuntimeError(
            f"GenBank directory was not provided and no genbank_dir was found for pipeline '{args.pipeline_name}'"
        )

    fasta_dir = Path(fasta_dir_value).expanduser()
    genbank_dir = Path(genbank_dir_value).expanduser()
    output_csv = (
        Path(args.output_csv).expanduser()
        if args.output_csv
        else script_temp_dir() / "Accessions_metadata.csv"
    )
    workers = max(1, args.workers)

    if not fasta_dir.is_dir():
        raise RuntimeError(f"FASTA directory was not found: {fasta_dir}")
    if not genbank_dir.is_dir():
        raise RuntimeError(f"GenBank directory was not found: {genbank_dir}")
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    fasta_paths = iter_fasta_paths(fasta_dir)
    if not fasta_paths:
        raise RuntimeError(f"No FASTA files were found in directory: {fasta_dir}")
    target_accessions, _fasta_file_by_accession, refid_by_accession, total_accession_records = collect_accessions_from_fastas(fasta_paths)

    extracted_dir = Path(tempfile.mkdtemp(prefix="accessions_metadata_", dir=script_temp_dir()))
    try:
        gb_archive_file_by_accession = extract_matching_genbank_records(
            genbank_dir,
            set(target_accessions),
            extracted_dir,
        )
        rows, missing_accessions = build_rows(
            extracted_dir,
            target_accessions,
            refid_by_accession,
            gb_archive_file_by_accession,
            workers,
        )
        rows.sort(key=lambda row: (row["RefID"], row["Accession"]))
        if len(rows) != total_accession_records:
            raise RuntimeError(
                f"Output row count ({len(rows)}) does not match total FASTA accession record count ({total_accession_records})"
            )
        fieldnames = ordered_fieldnames(rows)
        write_csv(output_csv, rows, fieldnames)

        print(
            {
                "output_csv": str(output_csv.resolve()),
                "pipeline_name": args.pipeline_name,
                "fasta_dir": str(fasta_dir.resolve()),
                "genbank_dir": str(genbank_dir.resolve()),
                "fasta_file_count": len(fasta_paths),
                "fasta_accession_record_count": total_accession_records,
                "accession_count": len(target_accessions),
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
