#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO


FASTA_EXTENSIONS = {".fa", ".faa", ".fasta", ".fna", ".fas", ".ffn", ".frn", ".seq"}
TARGET_GENES = ("NS3", "NS5A", "NS5B")
NUCLEOTIDE_CHARS = set("ACGTUNWSMKRYBDHV")
BLAST_OUTFMT = "6 qseqid sseqid bitscore evalue length pident"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read all FASTA records from a folder, detect whether each accession has NS3, NS5A, "
            "and NS5B by BLAST against HCV.fasta, and write a CSV summary."
        )
    )
    parser.add_argument(
        "--fasta-dir",
        required=True,
        help="Directory containing FASTA files",
    )
    parser.add_argument(
        "--hcv-fasta",
        default="HCV.fasta",
        help="Path to HCV.fasta containing NS3, NS5A, and NS5B references",
    )
    parser.add_argument(
        "--output-csv",
        default="Outputs/accession_gene_hits.csv",
        help="Output CSV path",
    )
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def resolve_repo_path(value: str) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return repo_root() / path


def iter_fasta_paths(fasta_dir: Path) -> list[Path]:
    return sorted(
        path for path in fasta_dir.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    )


def normalize_sequence(sequence: str) -> str:
    return re.sub(r"[^A-Z*]", "", sequence.upper())


def sequence_is_nucleotide(sequence: str) -> bool:
    letters = [char for char in sequence if char.isalpha()]
    if not letters:
        return False
    return all(char in NUCLEOTIDE_CHARS for char in letters)


def collect_sequences(fasta_paths: list[Path]) -> dict[str, str]:
    sequences: dict[str, str] = {}
    for fasta_path in fasta_paths:
        for record in SeqIO.parse(fasta_path, "fasta"):
            accession = str(record.id).strip()
            if not accession:
                continue
            if accession in sequences:
                raise RuntimeError(f"Duplicate accession across FASTA files or records: {accession}")
            normalized = normalize_sequence(str(record.seq))
            if normalized:
                sequences[accession] = normalized
    if not sequences:
        raise RuntimeError("No accession sequences were found in the FASTA folder")
    return sequences


def write_fasta(path: Path, entries: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in entries:
            handle.write(f">{header}\n{sequence}\n")


def run_command(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Required executable not found: {cmd[0]}") from exc
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() if exc.stderr else ""
        stdout = exc.stdout.strip() if exc.stdout else ""
        details = stderr or stdout or str(exc)
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\n{details}") from exc


def build_blast_db(reference_fasta: Path, db_prefix: Path) -> None:
    run_command(
        [
            "makeblastdb",
            "-in",
            str(reference_fasta),
            "-dbtype",
            "prot",
            "-out",
            str(db_prefix),
        ]
    )


def parse_blast_hits(blast_output: str) -> dict[str, set[str]]:
    hits_by_accession: dict[str, set[str]] = {}
    for line in blast_output.splitlines():
        if not line.strip():
            continue
        columns = line.split("\t")
        if len(columns) < 2:
            continue
        accession = columns[0].strip()
        gene = columns[1].strip()
        if gene not in TARGET_GENES:
            continue
        hits_by_accession.setdefault(accession, set()).add(gene)
    return hits_by_accession


def run_blast(query_fasta: Path, db_prefix: Path, program: str) -> dict[str, set[str]]:
    if query_fasta.stat().st_size == 0:
        return {}
    completed = run_command(
        [
            program,
            "-query",
            str(query_fasta),
            "-db",
            str(db_prefix),
            "-outfmt",
            BLAST_OUTFMT,
            "-evalue",
            "1e-8",
        ]
    )
    return parse_blast_hits(completed.stdout)


def merge_hits(*hit_maps: dict[str, set[str]]) -> dict[str, set[str]]:
    merged: dict[str, set[str]] = {}
    for hit_map in hit_maps:
        for accession, genes in hit_map.items():
            merged.setdefault(accession, set()).update(genes)
    return merged


def detect_gene_hits(sequences: dict[str, str], reference_fasta: Path, temp_dir: Path) -> dict[str, set[str]]:
    protein_entries: list[tuple[str, str]] = []
    nucleotide_entries: list[tuple[str, str]] = []

    for accession, sequence in sorted(sequences.items()):
        if sequence_is_nucleotide(sequence):
            nucleotide_entries.append((accession, sequence))
        else:
            protein_entries.append((accession, sequence))

    protein_query = temp_dir / "queries_protein.fasta"
    nucleotide_query = temp_dir / "queries_nucleotide.fasta"
    write_fasta(protein_query, protein_entries)
    write_fasta(nucleotide_query, nucleotide_entries)

    db_prefix = temp_dir / "hcv_refs"
    build_blast_db(reference_fasta, db_prefix)
    protein_hits = run_blast(protein_query, db_prefix, "blastp")
    nucleotide_hits = run_blast(nucleotide_query, db_prefix, "blastx")
    return merge_hits(protein_hits, nucleotide_hits)


def write_csv(path: Path, accessions: list[str], hits_by_accession: dict[str, set[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Accession", "NS3", "NS5A", "NS5B"])
        writer.writeheader()
        for accession in accessions:
            hits = hits_by_accession.get(accession, set())
            writer.writerow(
                {
                    "Accession": accession,
                    "NS3": "yes" if "NS3" in hits else "",
                    "NS5A": "yes" if "NS5A" in hits else "",
                    "NS5B": "yes" if "NS5B" in hits else "",
                }
            )


def main() -> int:
    args = parse_args()
    fasta_dir = resolve_repo_path(args.fasta_dir)
    reference_fasta = resolve_repo_path(args.hcv_fasta)
    output_csv = resolve_repo_path(args.output_csv)

    if not fasta_dir.is_dir():
        raise RuntimeError(f"FASTA directory not found: {fasta_dir}")
    if not reference_fasta.is_file():
        raise RuntimeError(f"HCV reference FASTA not found: {reference_fasta}")

    fasta_paths = iter_fasta_paths(fasta_dir)
    if not fasta_paths:
        raise RuntimeError(f"No FASTA files found in directory: {fasta_dir}")

    sequences = collect_sequences(fasta_paths)
    accessions = sorted(sequences)

    with tempfile.TemporaryDirectory(prefix="detect_accession_hcv_genes_") as temp_dir_text:
        hits_by_accession = detect_gene_hits(sequences, reference_fasta, Path(temp_dir_text))

    write_csv(output_csv, accessions, hits_by_accession)
    print(output_csv.relative_to(repo_root()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
