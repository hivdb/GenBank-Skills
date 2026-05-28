#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

try:
    import tomllib
except ModuleNotFoundError as exc:  # pragma: no cover
    raise SystemExit(f"tomllib is required: {exc}") from exc


FASTA_EXTENSIONS = {".fa", ".faa", ".fasta", ".fna", ".fas", ".ffn", ".frn", ".seq"}
GENBANK_EXTENSIONS = {".gb", ".gbk", ".genbank", ".seq"}
ACCESSION_RE = re.compile(r"^ACCESSION\s+(.+)$", re.MULTILINE)
VERSION_RE = re.compile(r"^VERSION\s+(\S+)", re.MULTILINE)


@dataclass(frozen=True)
class GenbankRecordLocation:
    source_path: Path
    byte_start: int
    byte_end: int


GenbankIndex = dict[str, GenbankRecordLocation]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "For each FASTA file in a directory, write one combined GenBank file "
            "containing matching records from a local GenBank directory."
        )
    )
    parser.add_argument(
        "--fasta-dir",
        "--fast-dir",
        default="",
        dest="fasta_dir",
        help="Directory containing FASTA files; defaults to [ns3].fasta_pool in pipeline.local.toml",
    )
    parser.add_argument(
        "--genbank-dir",
        default="",
        help="Directory containing local GenBank flatfiles; defaults to [ns3].genbank_dir in pipeline.local.toml",
    )
    parser.add_argument(
        "--pipeline-name",
        default="ns3",
        help="Pipeline section in pipeline.local.toml to use for default fasta/genbank paths",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs/genbank_by_fasta",
        help="Directory where combined GenBank files will be written",
    )
    parser.add_argument(
        "--rebuild-index",
        action="store_true",
        help="Re-scan the GenBank directory and rebuild the accession index even if it already exists",
    )
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def default_temp_dir(output_dir: Path) -> Path:
    return repo_root() / "temp" / output_dir.name


def convert_tsv_to_csv(source_path: Path, dest_path: Path) -> None:
    with source_path.open(encoding="utf-8", newline="") as source_handle:
        reader = csv.reader(source_handle, delimiter="\t")
        with dest_path.open("w", encoding="utf-8", newline="") as dest_handle:
            writer = csv.writer(dest_handle)
            writer.writerows(reader)


def migrate_legacy_tsv_indexes(output_dir: Path, temp_dir: Path) -> None:
    migrations = {
        output_dir / "accession_index.tsv": output_dir / "accession_index.csv",
        output_dir / "genbank_record_index.tsv": temp_dir / "genbank_record_index.csv",
        output_dir / "genbank_file_index.tsv": temp_dir / "genbank_file_index.csv",
        output_dir / "genbank_accession_source_index.tsv": temp_dir / "genbank_accession_source_index.csv",
        output_dir / "output_accession_index.tsv": temp_dir / "output_accession_index.csv",
    }
    for source_path, dest_path in migrations.items():
        if source_path.is_file():
            convert_tsv_to_csv(source_path, dest_path)
            source_path.unlink()


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

    defaults: dict[str, str] = {}
    fasta_pool = merged.get("fasta_pool")
    genbank_dir = merged.get("genbank_dir")
    if fasta_pool:
        defaults["fasta_dir"] = str(resolve_config_path(str(fasta_pool)))
    if genbank_dir:
        defaults["genbank_dir"] = str(resolve_config_path(str(genbank_dir)))
    return defaults


def iter_fasta_paths(fasta_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in fasta_dir.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    )


def collect_fasta_accessions(path: Path) -> list[str]:
    accessions: list[str] = []
    seen: set[str] = set()
    for record in SeqIO.parse(path, "fasta"):
        accession = str(record.id).strip()
        if not accession or accession in seen:
            continue
        accessions.append(accession)
        seen.add(accession)
    return accessions


def iter_genbank_paths(genbank_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in genbank_dir.iterdir()
        if path.is_file() and path.suffix.lower() in GENBANK_EXTENSIONS
    )


def iter_raw_genbank_records_with_offsets(path: Path):
    chunks: list[bytes] = []
    byte_start = 0
    with path.open("rb") as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            chunks.append(line)
            if line.rstrip(b"\r\n") == b"//":
                byte_end = handle.tell()
                yield byte_start, byte_end, b"".join(chunks).decode("utf-8", errors="replace")
                chunks = []
                byte_start = byte_end
    if chunks:
        byte_end = byte_start + sum(len(chunk) for chunk in chunks)
        tail = b"".join(chunks).strip()
        if tail:
            yield byte_start, byte_end, b"".join(chunks).decode("utf-8", errors="replace")


def accessions_from_raw_record(record_text: str) -> set[str]:
    accessions: set[str] = set()
    accession_match = ACCESSION_RE.search(record_text)
    if accession_match:
        accessions.update(accession_match.group(1).split())

    version_match = VERSION_RE.search(record_text)
    if version_match:
        version = version_match.group(1).strip()
        accessions.add(version)
        accessions.add(version.split(".", 1)[0])

    return accessions


def build_genbank_record_index_file(
    genbank_dir: Path,
    index_path: Path,
) -> dict[Path, set[str]]:
    accessions_by_file: dict[Path, set[str]] = {}
    seen_accessions: set[str] = set()
    with index_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["accession", "genbank_file", "byte_start", "byte_end", "is_primary_match"],
        )
        writer.writeheader()
        for genbank_path in iter_genbank_paths(genbank_dir):
            file_accessions: set[str] = set()
            for byte_start, byte_end, record_text in iter_raw_genbank_records_with_offsets(genbank_path):
                for accession in accessions_from_raw_record(record_text):
                    file_accessions.add(accession)
                    is_primary_match = "yes" if accession not in seen_accessions else "no"
                    seen_accessions.add(accession)
                    writer.writerow(
                        {
                            "accession": accession,
                            "genbank_file": str(genbank_path),
                            "byte_start": byte_start,
                            "byte_end": byte_end,
                            "is_primary_match": is_primary_match,
                        }
                    )
            accessions_by_file[genbank_path] = file_accessions
    return accessions_by_file


def load_genbank_record_index(index_path: Path) -> tuple[GenbankIndex, dict[Path, set[str]]]:
    index: GenbankIndex = {}
    accessions_by_file: dict[Path, set[str]] = {}
    with index_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            accession = row["accession"]
            source_path = Path(row["genbank_file"])
            byte_start = int(row["byte_start"])
            byte_end = int(row["byte_end"])
            accessions_by_file.setdefault(source_path, set()).add(accession)
            if accession not in index:
                index[accession] = GenbankRecordLocation(source_path, byte_start, byte_end)
    return index, accessions_by_file


def ensure_genbank_record_index(
    genbank_dir: Path,
    temp_dir: Path,
    rebuild: bool,
) -> tuple[Path, GenbankIndex, dict[Path, set[str]]]:
    index_path = temp_dir / "genbank_record_index.csv"
    if rebuild or not index_path.is_file():
        build_genbank_record_index_file(genbank_dir, index_path)
    genbank_index, accessions_by_file = load_genbank_record_index(index_path)
    return index_path, genbank_index, accessions_by_file


def write_genbank_source_indexes(
    temp_dir: Path,
    accessions_by_file: dict[Path, set[str]],
) -> None:
    file_index_path = temp_dir / "genbank_file_index.csv"
    with file_index_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["genbank_file", "accession_count", "accessions"],
        )
        writer.writeheader()
        for genbank_path in sorted(accessions_by_file):
            accessions = sorted(accessions_by_file[genbank_path])
            writer.writerow(
                {
                    "genbank_file": str(genbank_path),
                    "accession_count": len(accessions),
                    "accessions": ";".join(accessions),
                }
            )

    accession_index_path = temp_dir / "genbank_accession_source_index.csv"
    with accession_index_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["accession", "genbank_file"],
        )
        writer.writeheader()
        for genbank_path in sorted(accessions_by_file):
            for accession in sorted(accessions_by_file[genbank_path]):
                writer.writerow({"accession": accession, "genbank_file": str(genbank_path)})


def output_path_for_fasta(output_dir: Path, fasta_path: Path) -> Path:
    return output_dir / f"{fasta_path.stem}.gb"


def write_records_for_fasta(
    output_path: Path,
    accessions: list[str],
    genbank_index: GenbankIndex,
) -> tuple[int, list[dict[str, str]]]:
    found_records: list[bytes] = []
    index_rows: list[dict[str, str]] = []
    for accession in accessions:
        matched_accession = accession
        location = genbank_index.get(accession)
        if location is None and "." in accession:
            matched_accession = accession.split(".", 1)[0]
            location = genbank_index.get(matched_accession)
        if location is not None:
            with location.source_path.open("rb") as handle:
                handle.seek(location.byte_start)
                record_text = handle.read(location.byte_end - location.byte_start)
            found_records.append(record_text)
            index_rows.append(
                {
                    "accession": accession,
                    "matched_accession": matched_accession,
                    "output_genbank_file": str(output_path),
                    "source_genbank_file": str(location.source_path),
                    "byte_start": str(location.byte_start),
                    "byte_end": str(location.byte_end),
                    "found": "yes",
                }
            )
        else:
            index_rows.append(
                {
                    "accession": accession,
                    "matched_accession": "",
                    "output_genbank_file": str(output_path),
                    "source_genbank_file": "",
                    "byte_start": "",
                    "byte_end": "",
                    "found": "no",
                }
            )

    output_path.write_bytes(b"".join(found_records))
    return len(found_records), index_rows


def write_output_accession_index(output_dir: Path, rows: list[dict[str, str]]) -> None:
    path = output_dir / "accession_index.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["accession", "output_genbank_file"],
        )
        writer.writeheader()
        for row in rows:
            if row["found"] == "yes":
                writer.writerow(
                    {
                        "accession": row["accession"],
                        "output_genbank_file": Path(row["output_genbank_file"]).name,
                    }
                )


def write_detailed_output_accession_index(temp_dir: Path, rows: list[dict[str, str]]) -> None:
    path = temp_dir / "output_accession_index.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "fasta_file",
                "accession",
                "found",
                "matched_accession",
                "output_genbank_file",
                "source_genbank_file",
                "byte_start",
                "byte_end",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    **row,
                    "output_genbank_file": Path(row["output_genbank_file"]).name
                    if row["output_genbank_file"]
                    else "",
                    "source_genbank_file": Path(row["source_genbank_file"]).name
                    if row["source_genbank_file"]
                    else "",
                }
            )


def main() -> None:
    args = parse_args()
    defaults = load_pipeline_defaults(args.pipeline_name)
    fasta_dir_text = args.fasta_dir or defaults.get("fasta_dir", "")
    genbank_dir_text = args.genbank_dir or defaults.get("genbank_dir", "")
    if not fasta_dir_text:
        raise SystemExit(
            f"--fasta-dir was not provided and [{args.pipeline_name}].fasta_pool was not found in pipeline.local.toml"
        )
    if not genbank_dir_text:
        raise SystemExit(
            f"--genbank-dir was not provided and [{args.pipeline_name}].genbank_dir was not found in pipeline.local.toml"
        )

    fasta_dir = Path(fasta_dir_text).expanduser()
    genbank_dir = Path(genbank_dir_text).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    temp_dir = default_temp_dir(output_dir)

    if not fasta_dir.is_dir():
        raise SystemExit(f"FASTA directory was not found: {fasta_dir}")
    if not genbank_dir.is_dir():
        raise SystemExit(f"GenBank directory was not found: {genbank_dir}")

    fasta_paths = iter_fasta_paths(fasta_dir)
    if not fasta_paths:
        raise SystemExit(f"No FASTA files were found in: {fasta_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)
    migrate_legacy_tsv_indexes(output_dir, temp_dir)
    index_path, genbank_index, accessions_by_file = ensure_genbank_record_index(
        genbank_dir,
        temp_dir,
        args.rebuild_index,
    )
    write_genbank_source_indexes(temp_dir, accessions_by_file)

    print(f"genbank_record_index\t{index_path}")
    print("fasta_file\taccession_count\tfound_in_genbank_dir\toutput_file")
    output_index_rows: list[dict[str, str]] = []
    for fasta_path in fasta_paths:
        accessions = collect_fasta_accessions(fasta_path)
        output_path = output_path_for_fasta(output_dir, fasta_path)
        found_count, index_rows = write_records_for_fasta(output_path, accessions, genbank_index)
        for row in index_rows:
            row["fasta_file"] = fasta_path.name
        output_index_rows.extend(index_rows)
        print(f"{fasta_path.name}\t{len(accessions)}\t{found_count}\t{output_path}")
    write_output_accession_index(output_dir, output_index_rows)
    write_detailed_output_accession_index(temp_dir, output_index_rows)


if __name__ == "__main__":
    main()
