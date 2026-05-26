#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


TARGET_GENES = ("NS3", "NS5A_NTD", "NS5B")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract genotype-level HCV amino-acid references for NS3, NS5A_NTD, "
            "and NS5B from HCV_GT_Refs_By_Gene_AA.json and write them as FASTA."
        )
    )
    parser.add_argument(
        "--input-json",
        default="HCV_GT_Refs_By_Gene_AA.json",
        help="Path to HCV_GT_Refs_By_Gene_AA.json",
    )
    parser.add_argument(
        "--output-fasta",
        default="HCV_GT_Refs_NS3_NS5A_NTD_NS5B_AA.fasta",
        help="Path to output FASTA file",
    )
    parser.add_argument(
        "--ns5a-label",
        choices=("NS5A_NTD", "NS5A"),
        default="NS5A_NTD",
        help="Header label to use for NS5A_NTD entries in the FASTA output",
    )
    return parser.parse_args()


def wrap_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def load_entries(json_path: Path, ns5a_label: str) -> list[tuple[int, str, str]]:
    with json_path.open(encoding="utf-8") as handle:
        records = json.load(handle)

    entries: list[tuple[int, str, str]] = []
    seen: set[tuple[int, str]] = set()
    for record in records:
        gene = str(record.get("abstractGene", "")).strip()
        if gene not in TARGET_GENES:
            continue

        strain = str(record.get("strain", "")).strip()
        if not strain.startswith("HCV"):
            continue
        gt_text = strain.removeprefix("HCV")
        if not gt_text.isdigit():
            continue
        genotype = int(gt_text)

        header_gene = ns5a_label if gene == "NS5A_NTD" else gene
        key = (genotype, header_gene)
        if key in seen:
            continue

        sequence = str(record.get("refSequence", "")).strip().upper()
        if not sequence:
            continue

        seen.add(key)
        entries.append((genotype, header_gene, sequence))

    entries.sort(key=lambda item: (item[0], TARGET_GENES.index("NS5A_NTD" if item[1] in {"NS5A", "NS5A_NTD"} else item[1])))
    return entries


def main() -> None:
    args = parse_args()
    json_path = Path(args.input_json).expanduser().resolve()
    output_path = Path(args.output_fasta).expanduser().resolve()

    entries = load_entries(json_path, args.ns5a_label)
    if not entries:
        raise SystemExit(f"No matching entries found in {json_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for genotype, gene, sequence in entries:
            handle.write(f">genotype {genotype} | {gene}\n")
            handle.write(f"{wrap_sequence(sequence)}\n")
            print(
                f"genotype {genotype} | {gene} | length {len(sequence)} | begin {sequence[:3]} | end {sequence[-3:]}"
            )

    print(output_path)


if __name__ == "__main__":
    main()
