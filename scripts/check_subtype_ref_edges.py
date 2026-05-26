#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


DEFAULT_GT_FASTA = Path("HCV_GT_Refs_NS3_NS5A_NTD_NS5B_AA.fasta")
DEFAULT_SUBTYPE_FASTAS = (
    Path("temp/hcv_gene_subtype_refs/hcv_subtype_gene_refs_ns3_aa.fasta"),
    Path("temp/hcv_gene_subtype_refs/hcv_subtype_gene_refs_ns5a_ntd_aa.fasta"),
    Path("temp/hcv_gene_subtype_refs/hcv_subtype_gene_refs_ns5b_aa.fasta"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Check whether subtype AA references match the beginning and ending "
            "edge amino acids of the corresponding genotype AA references."
        )
    )
    parser.add_argument(
        "--gt-fasta",
        default=str(DEFAULT_GT_FASTA),
        help="Path to genotype AA FASTA created from HCV_GT_Refs_By_Gene_AA.json",
    )
    parser.add_argument(
        "--subtype-fasta",
        action="append",
        default=None,
        help=(
            "Subtype AA FASTA to check. Repeat to pass multiple files. "
            f"Defaults to: {', '.join(str(path) for path in DEFAULT_SUBTYPE_FASTAS)}"
        ),
    )
    parser.add_argument(
        "--edge-aa",
        type=int,
        choices=(2, 3),
        default=3,
        help="Number of amino acids to compare at the beginning and end",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    entries: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "".join(seq_parts).upper()))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        entries.append((header, "".join(seq_parts).upper()))
    return entries


def normalize_gene(gene: str) -> str:
    gene = gene.strip()
    return "NS5A_NTD" if gene == "NS5A" else gene


def parse_gt_header(header: str) -> tuple[str, str]:
    parts = [part.strip() for part in header.split("|")]
    if len(parts) != 2 or not parts[0].startswith("genotype "):
        raise ValueError(f"Unrecognized genotype FASTA header: {header}")
    genotype = parts[0].removeprefix("genotype ").strip()
    gene = normalize_gene(parts[1])
    return genotype, gene


def parse_subtype_header(header: str) -> dict[str, str]:
    fields: dict[str, str] = {}
    for part in header.split("|"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        fields[key.strip()] = value.strip()
    required = {"gene", "genotype", "subtype"}
    missing = sorted(required - fields.keys())
    if missing:
        raise ValueError(f"Subtype FASTA header missing {', '.join(missing)}: {header}")
    fields["gene"] = normalize_gene(fields["gene"])
    return fields


def load_gt_references(path: Path) -> dict[tuple[str, str], str]:
    refs: dict[tuple[str, str], str] = {}
    for header, sequence in read_fasta(path):
        key = parse_gt_header(header)
        refs[key] = sequence
    return refs


def compare_edges(
    subtype_path: Path,
    gt_refs: dict[tuple[str, str], str],
    edge_aa: int,
) -> int:
    mismatch_count = 0
    for header, subtype_seq in read_fasta(subtype_path):
        meta = parse_subtype_header(header)
        genotype = meta["genotype"]
        gene = meta["gene"]
        subtype = meta["subtype"]
        gt_key = (genotype, gene)
        if gt_key not in gt_refs:
            raise ValueError(
                f"Missing genotype reference for genotype={genotype}, gene={gene}"
            )
        gt_seq = gt_refs[gt_key]

        begin_gt = gt_seq[:edge_aa]
        begin_sub = subtype_seq[:edge_aa]
        end_gt = gt_seq[-edge_aa:]
        end_sub = subtype_seq[-edge_aa:]

        if begin_gt != begin_sub:
            mismatch_count += 1
            print(
                "\t".join(
                    [
                        "begin",
                        f"subtype={subtype}",
                        f"genotype={genotype}",
                        f"gene={gene}",
                        f"genotype_aa={begin_gt}",
                        f"subtype_aa={begin_sub}",
                        f"source={subtype_path.name}",
                    ]
                )
            )
        if end_gt != end_sub:
            mismatch_count += 1
            print(
                "\t".join(
                    [
                        "end",
                        f"subtype={subtype}",
                        f"genotype={genotype}",
                        f"gene={gene}",
                        f"genotype_aa={end_gt}",
                        f"subtype_aa={end_sub}",
                        f"source={subtype_path.name}",
                    ]
                )
            )
    return mismatch_count


def main() -> None:
    args = parse_args()
    gt_fasta = Path(args.gt_fasta).expanduser().resolve()
    subtype_fastas = (
        [Path(item).expanduser().resolve() for item in args.subtype_fasta]
        if args.subtype_fasta
        else [path.resolve() for path in DEFAULT_SUBTYPE_FASTAS]
    )

    gt_refs = load_gt_references(gt_fasta)
    mismatch_count = 0
    for subtype_path in subtype_fastas:
        mismatch_count += compare_edges(subtype_path, gt_refs, args.edge_aa)

    if mismatch_count:
        print(f"total_mismatches={mismatch_count}")


if __name__ == "__main__":
    main()
