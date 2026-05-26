#!/usr/bin/env python3
from __future__ import annotations

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
    rows: list[dict[str, str]] = []

    for gene, fasta_path in INPUT_FASTAS.items():
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

    with OUTPUT_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Gene", "Subtype", "beginAA", "endAA"])
        writer.writeheader()
        writer.writerows(rows)

    print(OUTPUT_CSV.relative_to(REPO_ROOT))


if __name__ == "__main__":
    main()
