#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split included_accessions_metadata.csv into one CSV per RefID, applying "
            "RefID-specific row filters and writing a filter summary."
        )
    )
    parser.add_argument("--input-csv", required=True, help="Path to included_accessions_metadata.csv")
    parser.add_argument("--output-dir", required=True, help="Directory for per-RefID CSV files")
    return parser.parse_args()


def sanitize_filename(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return text.strip("._-") or "unknown"


def load_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if "RefID" not in fieldnames:
            raise RuntimeError(f"Column 'RefID' was not found in {path}")
        return fieldnames, list(reader)


def load_accessions(path: Path) -> set[str]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        accession_column = "Accession" if "Accession" in fieldnames else "IsolateID"
        if accession_column not in fieldnames:
            raise RuntimeError(f"Column 'Accession' or 'IsolateID' was not found in {path}")
        return {
            accession
            for row in reader
            if (accession := (row.get(accession_column) or "").strip())
        }


def refid_filter_description(refid: str) -> str:
    if refid == "30":
        return "source_isolate contains Day1"
    if refid == "85":
        return "Accession in 85.csv"
    if refid == "142":
        return "source_isolate contains baseline"
    if refid == "192":
        return "source_isolate contains day 1"
    if refid == "346":
        return "source_isolate contains baseline/D0"
    if refid == "499":
        return "source_isolate contains HCC"
    if refid == "600":
        return "source_isolate does not contain failure"
    if refid == "661":
        return "source_isolation_source == plasma"
    if refid == "884":
        return "source_isolate contains Pre-TH"
    if refid == "943":
        return "source_isolate contains Day 1"
    if refid == "1356":
        return "source_isolate does not contain IC"
    if refid == "2008":
        return "source_isolate does not contain chimpanzee"
    if refid == "2110":
        return "source_isolate contains T0"
    if refid == "2116":
        return "source_collection_date before 2011"
    if refid == "2138":
        return "source_isolate contains Week 0"
    if refid == "2150":
        return "source_isolate contains b"
    if refid == "2168":
        return "source_isolate contains pre"
    if refid == "2178":
        return "source_isolation_source == plasma"
    if refid == "2227":
        return "Accession in 2227_Nguyen_(2015)_w_metadata_filtered.csv"
    raise KeyError(refid)


def text_contains(row: dict[str, str], column: str, needle: str) -> bool:
    return needle.casefold() in (row.get(column) or "").casefold()


def text_does_not_contain(row: dict[str, str], column: str, needle: str) -> bool:
    return needle.casefold() not in (row.get(column) or "").casefold()


def text_equals(row: dict[str, str], column: str, value: str) -> bool:
    return (row.get(column) or "").strip().casefold() == value.casefold()


def first_year(value: str) -> int | None:
    match = re.search(r"\b(19|20)\d{2}\b", value)
    return int(match.group(0)) if match else None


def row_is_kept(refid: str, row: dict[str, str]) -> bool:
    if refid == "30":
        return text_contains(row, "source_isolate", "Day1")
    if refid == "142":
        return text_contains(row, "source_isolate", "baseline")
    if refid == "192":
        return text_contains(row, "source_isolate", "day 1")
    if refid == "346":
        return text_contains(row, "source_isolate", "baseline/D0")
    if refid == "499":
        return text_contains(row, "source_isolate", "HCC")
    if refid == "600":
        return text_does_not_contain(row, "source_isolate", "failure")
    if refid == "661":
        return text_equals(row, "source_isolation_source", "plasma")
    if refid == "884":
        return text_contains(row, "source_isolate", "Pre-TH")
    if refid == "943":
        return text_contains(row, "source_isolate", "Day 1")
    if refid == "1356":
        return text_does_not_contain(row, "source_isolate", "IC")
    if refid == "2008":
        return text_does_not_contain(row, "source_isolate", "chimpanzee")
    if refid == "2110":
        return text_contains(row, "source_isolate", "T0")
    if refid == "2116":
        year = first_year(row.get("source_collection_date") or "")
        return year is not None and year < 2011
    if refid == "2138":
        return text_contains(row, "source_isolate", "Week 0")
    if refid == "2150":
        return text_contains(row, "source_isolate", "b")
    if refid == "2168":
        return text_contains(row, "source_isolate", "pre")
    if refid == "2178":
        return text_equals(row, "source_isolation_source", "plasma")
    raise KeyError(refid)


def filtered_refids() -> set[str]:
    return {"30", "85", "142", "192", "346", "499", "600", "661", "884", "943", "1356", "2008", "2110", "2116", "2138", "2150", "2168", "2178", "2227"}


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_csv = Path(args.input_csv).expanduser()
    output_dir = Path(args.output_dir).expanduser()

    if not input_csv.is_file():
        raise RuntimeError(f"Input CSV was not found: {input_csv}")
    output_dir.mkdir(parents=True, exist_ok=True)

    fieldnames, rows = load_rows(input_csv)
    accession_filters = {
        "85": load_accessions(input_csv.parent / "85.csv"),
        "2227": load_accessions(input_csv.parent / "2227_Nguyen_(2015)_w_metadata_filtered.csv")
    }
    rows_by_refid: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        refid = (row.get("RefID") or "").strip()
        if refid:
            rows_by_refid[refid].append(row)

    summary_rows: list[dict[str, str | int]] = []
    for refid in sorted(filtered_refids(), key=lambda value: (int(value) if value.isdigit() else 10**12, value)):
        ref_rows = rows_by_refid[refid]
        if refid in accession_filters:
            kept_rows = [
                row
                for row in ref_rows
                if (row.get("Accession") or "").strip() in accession_filters[refid]
            ]
        else:
            kept_rows = [row for row in ref_rows if row_is_kept(refid, row)]
        output_csv = output_dir / f"RefID_{sanitize_filename(refid)}_metadata.csv"
        write_csv(output_csv, fieldnames, kept_rows)
        summary_rows.append(
            {
                "RefID": refid,
                "Filter": refid_filter_description(refid),
                "TotalRows": len(ref_rows),
                "KeptRows": len(kept_rows),
                "RemovedRows": len(ref_rows) - len(kept_rows),
                "OutputCSV": str(output_csv.resolve()),
            }
        )

    print(f"refid_count={len(summary_rows)}")
    print(f"input_row_count={len(rows)}")
    print(f"filtered_refids={','.join(sorted(filtered_refids()))}")
    print(f"output_dir={output_dir.resolve()}")
    for row in summary_rows:
        print(
            "filter_result="
            f"RefID:{row['RefID']},"
            f"Filter:{row['Filter']},"
            f"TotalRows:{row['TotalRows']},"
            f"KeptRows:{row['KeptRows']},"
            f"RemovedRows:{row['RemovedRows']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
