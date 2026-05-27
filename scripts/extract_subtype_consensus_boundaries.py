#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "outputs"
INPUT_FASTAS = {
    "NS3": OUTPUT_DIR / "NS3_Subtype_Consensus.fasta",
    "NS5A": OUTPUT_DIR / "NS5A_Subtype_Consensus.fasta",
    "NS5B": OUTPUT_DIR / "NS5B_Subtype_Consensus.fasta",
}
OUTPUT_CSV = OUTPUT_DIR / "Subtype_Consensus_Boundaries.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract first and last three amino acids from subtype consensus FASTA files."
    )
    parser.add_argument(
        "--input-root",
        default=str(OUTPUT_DIR),
        help=(
            "Directory containing gene subdirectories with subtype consensus FASTA files, "
            "or the FASTA files directly. Defaults to outputs/."
        ),
    )
    parser.add_argument(
        "--output-csv",
        default=str(OUTPUT_CSV),
        help="Path for the output CSV. Defaults to outputs/Subtype_Consensus_Boundaries.csv.",
    )
    return parser.parse_args()


def resolve_input_fastas(input_root: Path) -> dict[str, Path]:
    fastas: dict[str, Path] = {}
    for gene in INPUT_FASTAS:
        nested_path = input_root / gene / f"{gene}_Subtype_Consensus.fasta"
        direct_path = input_root / f"{gene}_Subtype_Consensus.fasta"
        fastas[gene] = nested_path if nested_path.exists() else direct_path
    return fastas


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks).strip().upper()))
            header = line[1:].strip()
            chunks = []
            continue
        chunks.append(line.strip())

    if header is not None:
        records.append((header, "".join(chunks).strip().upper()))
    return records


def parse_subtype(header: str) -> str:
    if "_" not in header:
        raise RuntimeError(f"Could not parse subtype from FASTA header: {header}")
    return header.split("_", 1)[1]


def main() -> None:
    args = parse_args()
    input_fastas = resolve_input_fastas(Path(args.input_root).expanduser())
    output_csv = Path(args.output_csv).expanduser()
    rows: list[dict[str, str]] = []

    for gene, fasta_path in input_fastas.items():
        if not fasta_path.exists():
            raise RuntimeError(f"Missing input FASTA: {fasta_path}")
        for header, sequence in read_fasta(fasta_path):
            if not sequence:
                continue
            rows.append(
                {
                    "Gene": gene,
                    "Subtype": parse_subtype(header),
                    "beginAA": sequence[:3],
                    "endAA": sequence[-3:],
                }
            )

    rows.sort(key=lambda row: (row["Gene"], row["Subtype"]))

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Gene", "Subtype", "beginAA", "endAA"])
        writer.writeheader()
        writer.writerows(rows)

    print(output_csv)


if __name__ == "__main__":
    main()
