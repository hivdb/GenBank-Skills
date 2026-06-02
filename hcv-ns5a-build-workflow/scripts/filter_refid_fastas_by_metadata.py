#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter copied RefID FASTA files in place using per-RefID metadata CSV "
            "files produced by split_refid_metadata_csv.py."
        )
    )
    parser.add_argument("--metadata-dir", required=True, help="Directory containing RefID_*_metadata.csv files")
    parser.add_argument("--fasta-dir", required=True, help="Directory containing copied RefID FASTA files")
    return parser.parse_args()


def refid_from_metadata_path(path: Path) -> str | None:
    match = re.fullmatch(r"RefID_(.+)_metadata\.csv", path.name)
    return match.group(1) if match else None


def refid_from_fasta_path(path: Path) -> str | None:
    match = re.match(r"([^_]+)_", path.name)
    return match.group(1) if match else None


def load_metadata_accessions(path: Path) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if "Accession" not in fieldnames:
            raise RuntimeError(f"Column 'Accession' was not found in {path}")
        return {
            accession
            for row in reader
            if (accession := (row.get("Accession") or "").strip())
        }


def read_fasta(path: Path) -> list[tuple[str, list[str]]]:
    records: list[tuple[str, list[str]]] = []
    header: str | None = None
    sequence_lines: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, sequence_lines))
                header = line.rstrip("\n")
                sequence_lines = []
            else:
                sequence_lines.append(line.rstrip("\n"))
    if header is not None:
        records.append((header, sequence_lines))
    return records


def header_accession(header: str) -> str:
    return header[1:].strip().split()[0]


def write_fasta(path: Path, records: list[tuple[str, list[str]]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence_lines in records:
            handle.write(f"{header}\n")
            for sequence_line in sequence_lines:
                handle.write(f"{sequence_line}\n")


def main() -> int:
    args = parse_args()
    metadata_dir = Path(args.metadata_dir).expanduser()
    fasta_dir = Path(args.fasta_dir).expanduser()

    if not metadata_dir.is_dir():
        raise RuntimeError(f"Metadata directory was not found: {metadata_dir}")
    if not fasta_dir.is_dir():
        raise RuntimeError(f"FASTA directory was not found: {fasta_dir}")

    metadata_paths = sorted(metadata_dir.glob("RefID_*_metadata.csv"))
    fasta_by_refid = {
        refid: path
        for path in sorted(fasta_dir.glob("*.fasta"))
        if (refid := refid_from_fasta_path(path))
    }

    total_before = 0
    total_after = 0
    filtered_refids = 0
    filtered_fasta_refids = 0
    missing_fasta_count = 0

    for metadata_path in metadata_paths:
        refid = refid_from_metadata_path(metadata_path)
        if not refid:
            continue
        filtered_refids += 1
        fasta_path = fasta_by_refid.get(refid)
        if fasta_path is None:
            missing_fasta_count += 1
            print(f"filter_result=RefID:{refid},Status:missing_fasta,TotalRecords:0,KeptRecords:0")
            continue

        allowed_accessions = load_metadata_accessions(metadata_path)
        records = read_fasta(fasta_path)
        kept_records = [
            record
            for record in records
            if header_accession(record[0]) in allowed_accessions
        ]
        write_fasta(fasta_path, kept_records)

        total_before += len(records)
        total_after += len(kept_records)
        filtered_fasta_refids += 1
        print(
            "filter_result="
            f"RefID:{refid},"
            f"TotalRecords:{len(records)},"
            f"KeptRecords:{len(kept_records)},"
            f"RemovedRecords:{len(records) - len(kept_records)}"
        )

    print(f"filtered_refid_count={filtered_refids}")
    print(f"filtered_fasta_refid_count={filtered_fasta_refids}")
    print(f"missing_fasta_count={missing_fasta_count}")
    print(f"total_records_before={total_before}")
    print(f"total_records_after={total_after}")
    print(f"total_records_removed={total_before - total_after}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
