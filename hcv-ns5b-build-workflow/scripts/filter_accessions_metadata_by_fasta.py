#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FASTA_EXTENSIONS = {".fa", ".faa", ".fasta", ".fna", ".fas", ".ffn", ".frn", ".seq"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collect accessions from included RefID FASTA files, filter Accessions_metadata.csv "
            "to those accessions, and report FASTA accessions missing from the metadata CSV."
        )
    )
    parser.add_argument("--fasta-dir", required=True, help="Directory containing included RefID FASTA files")
    parser.add_argument("--metadata-csv", required=True, help="Path to Accessions_metadata.csv")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for filtered metadata CSV and missing-accession reports",
    )
    parser.add_argument("--accession-column", default="Accession", help="Metadata CSV accession column")
    return parser.parse_args()


def looks_like_fasta(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS


def accession_from_header(line: str) -> str:
    return line[1:].strip().split()[0]


def collect_fasta_accessions(fasta_dir: Path) -> tuple[list[str], dict[str, list[str]]]:
    accessions: list[str] = []
    files_by_accession: dict[str, list[str]] = {}
    seen: set[str] = set()

    for fasta_path in sorted(path for path in fasta_dir.iterdir() if looks_like_fasta(path)):
        with fasta_path.open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line.startswith(">"):
                    continue
                accession = accession_from_header(line)
                if not accession:
                    continue
                files_by_accession.setdefault(accession, []).append(str(fasta_path))
                if accession not in seen:
                    seen.add(accession)
                    accessions.append(accession)
    return accessions, files_by_accession


def read_metadata_rows(
    metadata_csv: Path,
    accession_column: str,
) -> tuple[list[str], list[dict[str, str]], dict[str, list[dict[str, str]]]]:
    with metadata_csv.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if accession_column not in fieldnames:
            raise RuntimeError(f"Column '{accession_column}' was not found in {metadata_csv}")
        rows = list(reader)

    rows_by_accession: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        accession = (row.get(accession_column) or "").strip()
        if accession:
            rows_by_accession.setdefault(accession, []).append(row)
    return fieldnames, rows, rows_by_accession


def write_lines(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    fasta_dir = Path(args.fasta_dir).expanduser()
    metadata_csv = Path(args.metadata_csv).expanduser()
    output_dir = Path(args.output_dir).expanduser()

    if not fasta_dir.is_dir():
        raise RuntimeError(f"FASTA directory was not found: {fasta_dir}")
    if not metadata_csv.is_file():
        raise RuntimeError(f"Metadata CSV was not found: {metadata_csv}")
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_accessions, _files_by_accession = collect_fasta_accessions(fasta_dir)
    fasta_accession_set = set(fasta_accessions)
    fieldnames, metadata_rows, rows_by_accession = read_metadata_rows(metadata_csv, args.accession_column)

    filtered_rows: list[dict[str, str]] = []
    for row in metadata_rows:
        accession = (row.get(args.accession_column) or "").strip()
        if accession in fasta_accession_set:
            filtered_rows.append(row)

    metadata_accessions = set(rows_by_accession)
    missing_accessions = [accession for accession in fasta_accessions if accession not in metadata_accessions]

    filtered_csv = output_dir / "included_accessions_metadata.csv"
    missing_txt = output_dir / "missing_accessions_from_metadata.txt"

    write_csv(filtered_csv, fieldnames, filtered_rows)
    write_lines(missing_txt, missing_accessions)

    print(f"fasta_accession_count={len(fasta_accessions)}")
    print(f"metadata_accession_count={len(metadata_accessions)}")
    print(f"filtered_row_count={len(filtered_rows)}")
    print(f"missing_accession_count={len(missing_accessions)}")
    print(f"filtered_csv={filtered_csv.resolve()}")
    print(f"missing_accessions_file={missing_txt.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
