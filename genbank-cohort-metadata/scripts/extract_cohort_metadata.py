#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

PRIORITY_PERSON_FIELDS = [
    "patient_id",
    "patient",
    "person_id",
    "person",
    "subject_id",
    "subject",
    "individual",
    "individual_id",
    "donor",
    "donor_id",
    "host_subject_id",
    "specimen_voucher",
    "isolate",
    "strain",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract cohort metadata and person-level summary from GenBank accessions or a GenBank file."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--accession", action="append", help="GenBank accession; repeat for multiple accessions")
    source.add_argument("--gb-file", help="Local GenBank file path")
    parser.add_argument("--output-dir", default="outputs", help="Base output directory")
    parser.add_argument("--email", default="", help="Optional contact email for NCBI E-utilities")
    parser.add_argument("--tool", default="genbank-cohort-metadata", help="Tool name for NCBI E-utilities")
    parser.add_argument("--fully-rerun", action="store_true", help="Ignore cached GenBank files and redownload")
    return parser.parse_args()


def fetch_text(url: str) -> str:
    try:
        with urlopen(url) as response:
            return response.read().decode("utf-8")
    except HTTPError as exc:
        raise RuntimeError(f"Request failed with HTTP {exc.code}: {url}") from exc
    except URLError as exc:
        raise RuntimeError(f"Request failed for {url}: {exc.reason}") from exc


def fetch_genbank_record(accession: str, email: str, tool: str) -> str:
    query = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gbwithparts",
        "retmode": "text",
        "tool": tool,
    }
    if email:
        query["email"] = email
    text = fetch_text(f"{EFETCH_URL}?{urlencode(query)}")
    if not text.strip():
        raise RuntimeError(f"Empty GenBank response for {accession}")
    return text


def source_feature(record) -> SeqFeature | None:
    for feature in record.features:
        if feature.type == "source":
            return feature
    return None


def flatten_value(value: Any) -> str:
    if isinstance(value, list):
        return " | ".join(str(v) for v in value)
    return str(value)


def accession_from_record(record) -> str:
    if record.annotations.get("accessions"):
        return record.annotations["accessions"][0]
    return record.id.split(".")[0]


def sanitize_label(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return text.strip("._-") or "job"


def make_job_id(seed: str) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    digest = hashlib.sha1(seed.encode("utf-8")).hexdigest()[:8]
    return f"{timestamp}_{digest}"


def build_rows(records: Iterable[Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in records:
        source = source_feature(record)
        qualifiers = source.qualifiers if source else {}
        row: dict[str, str] = {
            "accession": accession_from_record(record),
            "version": record.id,
            "submission_date": str(record.annotations.get("date", "")),
            "definition": str(record.description or ""),
            "organism": str(record.annotations.get("organism", "")),
            "isolate": flatten_value(qualifiers.get("isolate", "")),
        }
        for key, value in sorted(qualifiers.items()):
            row[f"source_{key}"] = flatten_value(value)
        rows.append(row)
    return rows


def ordered_fieldnames(rows: list[dict[str, str]]) -> list[str]:
    preferred = [
        "accession",
        "version",
        "submission_date",
        "definition",
        "organism",
        "isolate",
    ]
    discovered = sorted({key for row in rows for key in row.keys() if key not in preferred})
    return preferred + discovered


def strip_clone_suffix(value: str) -> str:
    text = value.strip()
    patterns = [
        r"(?i)(.+?)[-_ ]clone[-_ ]?[a-z0-9]+$",
        r"(?i)(.+?)[-_ ]qs[-_ ]?[a-z0-9]+$",
        r"(?i)(.+?)[-_ ]cl[-_ ]?[a-z0-9]+$",
    ]
    for pattern in patterns:
        match = re.match(pattern, text)
        if match:
            return match.group(1).strip(" _-")
    return text


def likely_clone(row: dict[str, str], person_value: str) -> bool:
    text_pool = " ".join(
        str(row.get(key, ""))
        for key in ("isolate", "source_clone", "source_note", "definition")
    )
    if re.search(r"(?i)\bclone\b|\bquasispecies\b|\bqs\b", text_pool):
        return True
    isolate = row.get("isolate", "")
    if person_value and isolate and person_value != isolate:
        return True
    return False


def choose_person_field(rows: list[dict[str, str]]) -> tuple[str | None, dict[str, Any]]:
    candidates: list[tuple[int, str, int, int]] = []
    diagnostics: dict[str, Any] = {"candidate_scores": []}

    for field in PRIORITY_PERSON_FIELDS:
        column = field if field in rows[0] else f"source_{field}"
        values = [row.get(column, "").strip() for row in rows if row.get(column, "").strip()]
        if not values:
            continue
        distinct = {strip_clone_suffix(value) for value in values if strip_clone_suffix(value)}
        multi_accession_signal = sum(1 for count in Counter(strip_clone_suffix(value) for value in values).values() if count > 1)
        name_score = 5 if any(token in column.lower() for token in ("patient", "person", "subject", "individual", "donor")) else 1
        score = name_score * 100 + len(distinct) + multi_accession_signal * 10
        candidates.append((score, column, len(distinct), multi_accession_signal))
        diagnostics["candidate_scores"].append(
            {
                "field": column,
                "score": score,
                "distinct_ids": len(distinct),
                "multi_accession_signal": multi_accession_signal,
            }
        )

    if not candidates:
        return None, diagnostics

    candidates.sort(reverse=True)
    best = candidates[0][1]
    diagnostics["selected_field"] = best
    return best, diagnostics


def person_summary(rows: list[dict[str, str]], person_field: str | None) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if not person_field:
        return [], {"person_count": 0, "quasispecies_clone_count": 0}

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        raw = row.get(person_field, "").strip()
        if not raw:
            continue
        inferred = strip_clone_suffix(raw)
        if inferred:
            grouped[inferred].append(row)

    summary_rows: list[dict[str, Any]] = []
    clone_count = 0
    for person_id, person_rows in sorted(grouped.items()):
        accessions = [row["accession"] for row in person_rows]
        clone_flags = [likely_clone(row, person_id) for row in person_rows]
        has_clones = any(clone_flags) and len(accessions) > 1
        if has_clones:
            clone_count += 1
        summary_rows.append(
            {
                "person_id": person_id,
                "accession_count": len(accessions),
                "accessions": " | ".join(accessions),
                "likely_quasispecies_clones": "yes" if has_clones else "no",
            }
        )

    return summary_rows, {
        "person_count": len(summary_rows),
        "quasispecies_clone_count": clone_count,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def cohort_name_from_accessions(accessions: list[str]) -> str:
    cleaned = [acc.replace("/", "_") for acc in accessions]
    job_id = make_job_id("|".join(cleaned))
    if len(cleaned) <= 3:
        label = "_".join(cleaned)
    else:
        label = f"{cleaned[0]}_to_{cleaned[-1]}_{len(cleaned)}records"
    return f"cohort_{sanitize_label(label)}_{job_id}"


def cohort_name_from_gb_file(gb_file: Path) -> str:
    label = sanitize_label(gb_file.stem)
    job_id = make_job_id(str(gb_file.resolve()))
    return f"cohort_{label}_{job_id}"


def load_records_from_accessions(
    accessions: list[str],
    cohort_dir: Path,
    email: str,
    tool: str,
    fully_rerun: bool,
) -> list[Any]:
    records_dir = cohort_dir / "records"
    records_dir.mkdir(parents=True, exist_ok=True)
    records: list[Any] = []
    for accession in accessions:
        gb_path = records_dir / f"{accession.replace('/', '_')}.gb"
        if gb_path.exists() and not fully_rerun:
            text = gb_path.read_text(encoding="utf-8")
        else:
            text = fetch_genbank_record(accession, email, tool)
            gb_path.write_text(text, encoding="utf-8")
        records.extend(list(SeqIO.parse(gb_path, "genbank")))
    return records


def load_records_from_gb_file(gb_file: Path) -> list[Any]:
    if gb_file.is_dir():
        records: list[Any] = []
        for path in sorted(gb_file.glob("*.gb*")):
            records.extend(list(SeqIO.parse(path, "genbank")))
        return records
    return list(SeqIO.parse(gb_file, "genbank"))


def write_report(path: Path, cohort_name: str, row_count: int, person_field: str | None, person_meta: dict[str, Any]) -> None:
    lines = [
        "GenBank Cohort Metadata",
        "",
        f"Cohort: {cohort_name}",
        f"Record count: {row_count}",
        f"Selected person field: {person_field}",
        f"Number of persons: {person_meta.get('person_count')}",
        f"Likely quasispecies clone groups: {person_meta.get('quasispecies_clone_count')}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    base_output_dir = Path(args.output_dir)
    base_output_dir.mkdir(parents=True, exist_ok=True)
    job_id = ""
    gb_file: Path | None = None

    if args.accession:
        cohort_name = cohort_name_from_accessions(args.accession)
        job_id = "_".join(cohort_name.rsplit("_", 2)[-2:])
        cohort_dir = base_output_dir / cohort_name
        cohort_dir.mkdir(parents=True, exist_ok=True)
        records = load_records_from_accessions(args.accession, cohort_dir, args.email, args.tool, args.fully_rerun)
    else:
        gb_file = Path(args.gb_file).expanduser()
        cohort_name = cohort_name_from_gb_file(gb_file)
        job_id = "_".join(cohort_name.rsplit("_", 2)[-2:])
        cohort_dir = base_output_dir / cohort_name
        cohort_dir.mkdir(parents=True, exist_ok=True)
        records = load_records_from_gb_file(gb_file)

    if not records:
        raise RuntimeError("No GenBank records were found to process")

    rows = build_rows(records)
    metadata_path = cohort_dir / "metadata.csv"
    write_csv(metadata_path, rows, ordered_fieldnames(rows))

    person_field, diagnostics = choose_person_field(rows)
    summary_rows, person_meta = person_summary(rows, person_field)
    person_csv_path = cohort_dir / "person_summary.csv"
    write_csv(
        person_csv_path,
        summary_rows,
        ["person_id", "accession_count", "accessions", "likely_quasispecies_clones"],
    )

    summary = {
        "cohort_name": cohort_name,
        "record_count": len(rows),
        "selected_person_field": person_field,
        **person_meta,
        "metadata_csv": str(metadata_path),
        "person_summary_csv": str(person_csv_path),
        "input_mode": "accessions" if args.accession else "gb_file",
        "job_id": job_id,
        "source_accessions": args.accession or [],
        "source_gb_file": str(gb_file.resolve()) if not args.accession else None,
        "person_field_diagnostics": diagnostics,
    }
    write_json(cohort_dir / "summary.json", summary)
    write_report(cohort_dir / "summary_report.txt", cohort_name, len(rows), person_field, person_meta)

    print(json.dumps(summary, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
