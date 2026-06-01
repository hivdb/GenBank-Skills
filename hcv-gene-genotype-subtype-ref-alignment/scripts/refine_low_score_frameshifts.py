#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

NEG_INF = -10**9
MATCH_SCORE = 3
MISMATCH_SCORE = -2
STOP_SCORE = -6
AMBIG_SCORE = -3
FS1_PENALTY = -4
FS2_PENALTY = -5


@dataclass
class RefinedResult:
    alignment_score: float
    raw_dp_score: int
    frame: int
    nt_start_in_source: int
    nt_end_in_source: int
    query_aa_start: int
    query_aa_end: int
    refined_nt: str
    refined_aa: str
    match_line: str
    frameshift_events: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Refine low-scoring HCV gene subtype extractions with a local "
            "frameshift-aware codon alignment."
        )
    )
    parser.add_argument("--output-dir", required=True, help="Existing output directory from build_hcv_gene_subtype_refs.py")
    parser.add_argument("--gt-gene-aa-json", required=True, help="Path to HCV_GT_Refs_By_Gene_AA.json")
    parser.add_argument("--subtype-genome-na-json", required=True, help="Path to HCV_Subtype_Refs_By_Genome_NA.json")
    parser.add_argument("--threshold", type=float, default=0.80, help="Refine only records below this aa match proportion")
    parser.add_argument("--flank-nt", type=int, default=15, help="Extra nucleotides to include on each side of the initial window")
    return parser.parse_args()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def normalize_nt(sequence: str) -> str:
    return "".join(base for base in sequence.upper().replace("U", "T") if base in {"A", "C", "G", "T", "N"})


def parse_gt_from_name(name: str) -> str:
    prefix = "HCV"
    if not name.startswith(prefix):
        raise RuntimeError(f"Could not parse genotype from name: {name}")
    idx = len(prefix)
    while idx < len(name) and name[idx].isdigit():
        idx += 1
    return name[len(prefix):idx]


def build_gt_refs(rows: list[dict[str, Any]]) -> dict[tuple[str, str], str]:
    refs: dict[tuple[str, str], str] = {}
    for row in rows:
        gene = str(row["abstractGene"])
        if gene not in {"NS3", "NS5A_NTD", "NS5B"}:
            continue
        refs[(parse_gt_from_name(str(row["name"])), gene)] = str(row["refSequence"]).strip().upper()
    return refs


def translate_codon(codon: str) -> str:
    return CODON_TABLE.get(codon, "X")


def aa_pair_score(ref_aa: str, query_aa: str) -> int:
    if query_aa == "X":
        return AMBIG_SCORE
    if query_aa == "*":
        return STOP_SCORE
    return MATCH_SCORE if ref_aa == query_aa else MISMATCH_SCORE


def wrap(text: str, width: int = 80) -> list[str]:
    return [text[i : i + width] for i in range(0, len(text), width)]


def frameshift_refine(reference_aa: str, nt_window: str, source_nt_start: int) -> RefinedResult:
    m = len(reference_aa)
    n = len(nt_window)
    dp = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    trace: list[list[tuple[str, int, int, str] | None]] = [[None] * (n + 1) for _ in range(m + 1)]
    dp[0] = [0] * (n + 1)

    best_score = NEG_INF
    best_end_j = 0

    for i in range(0, m):
        ref_aa = reference_aa[i]
        for j in range(0, n + 1):
            current = dp[i][j]
            if current <= NEG_INF // 2:
                continue

            if j + 3 <= n:
                codon = nt_window[j : j + 3]
                aa = translate_codon(codon)
                score = current + aa_pair_score(ref_aa, aa)
                if score > dp[i + 1][j + 3]:
                    dp[i + 1][j + 3] = score
                    trace[i + 1][j + 3] = ("codon", i, j, aa)

            if j + 1 <= n:
                score = current + FS1_PENALTY
                if score > dp[i][j + 1]:
                    dp[i][j + 1] = score
                    trace[i][j + 1] = ("skip1", i, j, "")

            if j + 2 <= n:
                score = current + FS2_PENALTY
                if score > dp[i][j + 2]:
                    dp[i][j + 2] = score
                    trace[i][j + 2] = ("skip2", i, j, "")

    for j in range(0, n + 1):
        if dp[m][j] > best_score:
            best_score = dp[m][j]
            best_end_j = j

    if best_score <= NEG_INF // 2:
        raise RuntimeError("Frameshift refinement could not align the target window")

    aa_chars: list[str] = []
    match_chars: list[str] = []
    frameshift_events = 0
    kept_nt_chunks: list[str] = []
    i = m
    j = best_end_j
    while i > 0 or j > 0:
        step = trace[i][j]
        if step is None:
            break
        op, prev_i, prev_j, aa = step
        if op == "codon":
            aa_chars.append(aa)
            match_chars.append("|" if aa == reference_aa[i - 1] else " ")
            kept_nt_chunks.append(nt_window[prev_j:j])
        else:
            frameshift_events += 1
        i, j = prev_i, prev_j

    aa_chars.reverse()
    match_chars.reverse()
    kept_nt_chunks.reverse()
    refined_aa = "".join(aa_chars)
    refined_nt = "".join(kept_nt_chunks)
    matches = sum(1 for char in match_chars if char == "|")
    proportion = matches / len(reference_aa) if reference_aa else 0.0

    return RefinedResult(
        alignment_score=proportion,
        raw_dp_score=best_score,
        frame=-1,
        nt_start_in_source=source_nt_start + j + 1,
        nt_end_in_source=source_nt_start + best_end_j,
        query_aa_start=1,
        query_aa_end=len(reference_aa),
        refined_nt=refined_nt,
        refined_aa=refined_aa,
        match_line="".join(match_chars),
        frameshift_events=frameshift_events,
    )


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir).expanduser()
    gt_json = Path(args.gt_gene_aa_json).expanduser()
    subtype_json = Path(args.subtype_genome_na_json).expanduser()

    summary = load_json(output_dir / "summary.json")
    gt_rows = load_json(gt_json)
    subtype_rows = load_json(subtype_json)
    gt_refs = build_gt_refs(gt_rows)
    subtype_by_accession = {str(row["accession"]): row for row in subtype_rows}

    low_rows = [row for row in summary["records"] if float(row["alignment_score"]) < args.threshold]
    refined_rows: list[dict[str, Any]] = []

    csv_path = output_dir / "frameshift_refined_records_below_0_80.csv"
    md_by_gene = {
        "NS3": output_dir / "ns3_frameshift_refined_comparisons_below_0_80.txt",
        "NS5A_NTD": output_dir / "ns5a_ntd_frameshift_refined_comparisons_below_0_80.txt",
        "NS5B": output_dir / "ns5b_frameshift_refined_comparisons_below_0_80.txt",
    }
    md_lines = {
        gene: [
            f"{gene} Frameshift-Refined Comparisons",
            "",
            f"Threshold: aa_match_proportion < {args.threshold:.2f}",
            "",
        ]
        for gene in md_by_gene
    }

    for row in low_rows:
        accession = str(row["accession"])
        gene = str(row["gene"])
        genotype = str(row["genotype"])
        subtype_record = subtype_by_accession[accession]
        full_nt = normalize_nt(str(subtype_record["sequence"]))
        ref_aa = gt_refs[(genotype, gene)]

        initial_frame = int(row["frame"])
        initial_aa_start = int(row["query_aa_start"]) - 1
        initial_aa_end = int(row["query_aa_end"])
        initial_nt_start = initial_frame + initial_aa_start * 3
        initial_nt_end = initial_frame + initial_aa_end * 3
        window_start = max(0, initial_nt_start - args.flank_nt)
        window_end = min(len(full_nt), initial_nt_end + args.flank_nt)
        nt_window = full_nt[window_start:window_end]

        refined = frameshift_refine(ref_aa, nt_window, source_nt_start=window_start)
        refined_rows.append(
            {
                "gene": gene,
                "genotype": genotype,
                "subtype": row["subtype"],
                "genotype_name": row["genotype_name"],
                "accession": accession,
                "author_year": row["author_year"],
                "original_aa_match_proportion": row["alignment_score"],
                "refined_aa_match_proportion": refined.alignment_score,
                "improvement": refined.alignment_score - float(row["alignment_score"]),
                "refined_nt_start_in_source": refined.nt_start_in_source,
                "refined_nt_end_in_source": refined.nt_end_in_source,
                "frameshift_events": refined.frameshift_events,
                "raw_dp_score": refined.raw_dp_score,
            }
        )

        lines = md_lines[gene]
        lines.append(
            f"> comparison gene={gene} genotype={genotype} subtype={row['subtype']} accession={accession} "
            f"original={float(row['alignment_score']):.3f} refined={refined.alignment_score:.3f} "
            f"frameshifts={refined.frameshift_events}"
        )
        lines.append(
            f"source_nt_range: {refined.nt_start_in_source}-{refined.nt_end_in_source}"
        )
        lines.append("")
        for ref_chunk, match_chunk, tgt_chunk in zip(
            wrap(ref_aa),
            wrap(refined.match_line),
            wrap(refined.refined_aa),
        ):
            lines.append(f"REF  {ref_chunk}")
            lines.append(f"     {match_chunk}")
            lines.append(f"TGT  {tgt_chunk}")
            lines.append("")
        lines.append("-" * 100)
        lines.append("")

    with csv_path.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "gene",
            "genotype",
            "subtype",
            "genotype_name",
            "accession",
            "author_year",
            "original_aa_match_proportion",
            "refined_aa_match_proportion",
            "improvement",
            "refined_nt_start_in_source",
            "refined_nt_end_in_source",
            "frameshift_events",
            "raw_dp_score",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in refined_rows:
            writer.writerow(row)

    for gene, path in md_by_gene.items():
        path.write_text("\n".join(md_lines[gene]), encoding="utf-8")

    summary_payload = {
        "threshold": args.threshold,
        "refined_record_count": len(refined_rows),
        "csv": str(csv_path.resolve()),
        "comparison_files": {gene: str(path.resolve()) for gene, path in md_by_gene.items()},
    }
    (output_dir / "frameshift_refinement_summary.json").write_text(
        json.dumps(summary_payload, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary_payload, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
