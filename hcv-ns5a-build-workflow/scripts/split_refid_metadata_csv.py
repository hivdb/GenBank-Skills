#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split included_accessions_metadata.csv into per-RefID CSV files for "
            "NS5A-specific filters."
        )
    )
    parser.add_argument("--input-csv", required=True, help="Path to included_accessions_metadata.csv")
    parser.add_argument("--output-dir", required=True, help="Directory for per-RefID CSV files")
    return parser.parse_args()


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


def text_equals(row: dict[str, str], column: str, value: str) -> bool:
    return (row.get(column) or "").strip().casefold() == value.casefold()


def text_contains(row: dict[str, str], column: str, needle: str) -> bool:
    return needle.casefold() in (row.get(column) or "").casefold()


def text_does_not_contain(row: dict[str, str], column: str, needle: str) -> bool:
    return needle.casefold() not in (row.get(column) or "").casefold()


def refid_filter_description(refid: str) -> str:
    descriptions = {
        "17": "Accession in 17.csv",
        "29": "source_isolate contains SCRN",
        "50": "source_isolate contains week 0",
        "85": "Accession in 85.csv",
        "123": "source_isolate does not contain TF",
        "142": "source_isolate contains baseline",
        "165": "Accession in 165.csv",
        "192": "source_isolate contains day1",
        "288": "source_isolate contains pre",
        "346": "source_isolate contains baseline/D0",
        "535": "Accession in 535.csv",
        "600": "source_isolate does not contain failure",
        "661": "source_isolation_source == plasma",
    }
    return descriptions[refid]


def row_is_kept(refid: str, row: dict[str, str]) -> bool:
    if refid == "29":
        return text_contains(row, "source_isolate", "SCRN")
    if refid == "50":
        return text_contains(row, "source_isolate", "week 0")
    if refid == "123":
        return text_does_not_contain(row, "source_isolate", "TF")
    if refid == "142":
        return text_contains(row, "source_isolate", "baseline")
    if refid == "192":
        return text_contains(row, "source_isolate", "day1")
    if refid == "288":
        return text_contains(row, "source_isolate", "pre")
    if refid == "346":
        return text_contains(row, "source_isolate", "baseline/D0")
    if refid == "600":
        return text_does_not_contain(row, "source_isolate", "failure")
    if refid == "661":
        return text_equals(row, "source_isolation_source", "plasma")
    raise KeyError(refid)


def filtered_refids() -> set[str]:
    return {"17", "29", "50", "85", "123", "142", "165", "192", "288", "346", "535", "600", "661"}


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
        "17": load_accessions(input_csv.parent / "17.csv"),
        "85": load_accessions(input_csv.parent / "85.csv"),
        "165": load_accessions(input_csv.parent / "165.csv"),
        "535": load_accessions(input_csv.parent / "535.csv"),
    }
    rows_by_refid: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        refid = (row.get("RefID") or "").strip()
        if refid:
            rows_by_refid[refid].append(row)

    summary_rows: list[dict[str, str | int]] = []
    for refid in sorted(filtered_refids(), key=lambda value: (int(value), value)):
        ref_rows = rows_by_refid[refid]
        if refid in accession_filters:
            kept_rows = [
                row
                for row in ref_rows
                if (row.get("Accession") or "").strip() in accession_filters[refid]
            ]
        else:
            kept_rows = [row for row in ref_rows if row_is_kept(refid, row)]
        output_csv = output_dir / f"RefID_{refid}_metadata.csv"
        write_csv(output_csv, fieldnames, kept_rows)
        summary_rows.append(
            {
                "RefID": refid,
                "Filter": refid_filter_description(refid),
                "TotalRows": len(ref_rows),
                "KeptRows": len(kept_rows),
                "RemovedRows": len(ref_rows) - len(kept_rows),
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
