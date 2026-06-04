#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


FASTA_EXTENSIONS = {".fa", ".faa", ".fasta", ".fna", ".fas", ".ffn", ".frn", ".seq"}
GENOTYPE_SUBTYPE_COLUMNS = ["source_note", "source_organism", "source_strain", "source_serotype"]
GENOTYPE_TOKEN_RE = re.compile(r"\bgenotype\s*[:=]?\s*([1-8](?:[A-Za-z][A-Za-z0-9]*)?)\b", re.IGNORECASE)
SUBTYPE_TOKEN_RE = re.compile(r"\bsubtype\s*[:=]?\s*([1-8]?[A-Za-z][A-Za-z0-9]*)\b", re.IGNORECASE)
HCV_SUBTYPE_RE = re.compile(r"\bHCV[-\s]*([1-8][A-Za-z][A-Za-z0-9]*)\b", re.IGNORECASE)
BARE_GT_SUBTYPE_RE = re.compile(r"^\s*([1-8][A-Za-z][A-Za-z0-9]*)\s*$", re.IGNORECASE)


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


def genotype_from_token(token: str) -> str:
    match = re.match(r"([1-8])", token.strip())
    return match.group(1) if match else ""


def normalize_subtype(token: str, genotype: str) -> str:
    text = token.strip().lower()
    if re.fullmatch(r"[a-z][a-z0-9]*", text) and genotype:
        return f"{genotype}{text}"
    return text


def extract_genotype_subtype(value: str) -> tuple[str, str] | None:
    text = value.strip()
    if not text:
        return None

    genotype = ""
    subtype = ""

    genotype_match = GENOTYPE_TOKEN_RE.search(text)
    if genotype_match:
        token = genotype_match.group(1)
        genotype = genotype_from_token(token)
        if re.fullmatch(r"[1-8][A-Za-z][A-Za-z0-9]*", token):
            subtype = token.lower()

    subtype_match = SUBTYPE_TOKEN_RE.search(text)
    if subtype_match:
        subtype = normalize_subtype(subtype_match.group(1), genotype)
        genotype = genotype or genotype_from_token(subtype)

    hcv_match = HCV_SUBTYPE_RE.search(text)
    if hcv_match and not subtype:
        subtype = hcv_match.group(1).lower()
        genotype = genotype or genotype_from_token(subtype)

    bare_match = BARE_GT_SUBTYPE_RE.fullmatch(text)
    if bare_match and not (genotype or subtype):
        token = bare_match.group(1)
        genotype = genotype_from_token(token)
        if re.fullmatch(r"[1-8][A-Za-z][A-Za-z0-9]*", token):
            subtype = token.lower()

    if genotype or subtype:
        return genotype, subtype
    return None


def build_genotype_subtype_rows(
    metadata_rows: list[dict[str, str]],
    accession_column: str,
) -> list[dict[str, str]]:
    summary_rows: list[dict[str, str]] = []
    for row in metadata_rows:
        accession = (row.get(accession_column) or "").strip()
        summary = {
            "accession": accession,
            "genotype": "",
            "subtype": "",
            "column_name": "",
        }
        for column_name in GENOTYPE_SUBTYPE_COLUMNS:
            parsed = extract_genotype_subtype(row.get(column_name) or "")
            if parsed:
                summary["genotype"], summary["subtype"] = parsed
                summary["column_name"] = column_name
                break
        summary_rows.append(summary)
    return summary_rows


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
    genotype_subtype_csv = output_dir / "included_accessions_genotype_subtype.csv"
    missing_txt = output_dir / "missing_accessions_from_metadata.txt"

    write_csv(filtered_csv, fieldnames, filtered_rows)
    genotype_subtype_rows = build_genotype_subtype_rows(filtered_rows, args.accession_column)
    write_csv(genotype_subtype_csv, ["accession", "genotype", "subtype", "column_name"], genotype_subtype_rows)
    write_lines(missing_txt, missing_accessions)

    print(f"fasta_accession_count={len(fasta_accessions)}")
    print(f"metadata_accession_count={len(metadata_accessions)}")
    print(f"filtered_row_count={len(filtered_rows)}")
    print(f"missing_accession_count={len(missing_accessions)}")
    print(f"filtered_csv={filtered_csv.resolve()}")
    print(f"genotype_subtype_csv={genotype_subtype_csv.resolve()}")
    print(f"missing_accessions_file={missing_txt.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
