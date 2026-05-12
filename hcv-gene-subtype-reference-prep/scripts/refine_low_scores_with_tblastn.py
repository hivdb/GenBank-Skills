#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any


OUTFMT = "6 qseqid sseqid qlen qstart qend sstart send length nident mismatch gaps pident evalue bitscore qseq sseq frames"


@dataclass
class Hsp:
    qseqid: str
    sseqid: str
    qlen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    length: int
    nident: int
    mismatch: int
    gaps: int
    pident: float
    evalue: float
    bitscore: float
    qseq: str
    sseq: str
    frames: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Refine low-scoring HCV subtype references with tblastn segment chaining."
    )
    parser.add_argument("--output-dir", required=True, help="Output directory from build_hcv_gene_subtype_refs.py")
    parser.add_argument("--gt-gene-aa-json", required=True, help="Path to HCV_GT_Refs_By_Gene_AA.json")
    parser.add_argument("--subtype-genome-na-json", required=True, help="Path to HCV_Subtype_Refs_By_Genome_NA.json")
    parser.add_argument("--threshold", type=float, default=0.80, help="Only refine records below this score")
    return parser.parse_args()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def normalize_nt(sequence: str) -> str:
    return "".join(base for base in sequence.upper().replace("U", "T") if base in {"A", "C", "G", "T", "N"})


def parse_gt_from_name(name: str) -> str:
    if not name.startswith("HCV"):
        raise RuntimeError(f"Could not parse genotype from name: {name}")
    idx = 3
    while idx < len(name) and name[idx].isdigit():
        idx += 1
    return name[3:idx]


def build_gt_refs(rows: list[dict[str, Any]]) -> dict[tuple[str, str], str]:
    refs: dict[tuple[str, str], str] = {}
    for row in rows:
        gene = str(row["abstractGene"])
        if gene not in {"NS3", "NS5A_NTD", "NS5B"}:
            continue
        refs[(parse_gt_from_name(str(row["name"])), gene)] = str(row["refSequence"]).strip().upper()
    return refs


def write_fasta(path: Path, header: str, sequence: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{header}\n")
        for start in range(0, len(sequence), 70):
            handle.write(sequence[start : start + 70] + "\n")


def parse_hsp(line: str) -> Hsp:
    parts = line.rstrip("\n").split("\t")
    return Hsp(
        qseqid=parts[0],
        sseqid=parts[1],
        qlen=int(parts[2]),
        qstart=int(parts[3]),
        qend=int(parts[4]),
        sstart=int(parts[5]),
        send=int(parts[6]),
        length=int(parts[7]),
        nident=int(parts[8]),
        mismatch=int(parts[9]),
        gaps=int(parts[10]),
        pident=float(parts[11]),
        evalue=float(parts[12]),
        bitscore=float(parts[13]),
        qseq=parts[14],
        sseq=parts[15],
        frames=parts[16],
    )


def run_tblastn(query_aa: str, subject_nt: str, tag: str) -> list[Hsp]:
    with tempfile.TemporaryDirectory(prefix=f"tblastn_{tag}_") as tmpdir:
        tmp = Path(tmpdir)
        query_path = tmp / "query.faa"
        subject_path = tmp / "subject.fna"
        out_path = tmp / "hits.tsv"
        write_fasta(query_path, tag, query_aa)
        write_fasta(subject_path, tag, subject_nt)
        subprocess.run(
            [
                "tblastn",
                "-query",
                str(query_path),
                "-subject",
                str(subject_path),
                "-seg",
                "no",
                "-comp_based_stats",
                "F",
                "-max_hsps",
                "50",
                "-evalue",
                "1e-3",
                "-outfmt",
                OUTFMT,
                "-out",
                str(out_path),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        return [parse_hsp(line) for line in out_path.read_text(encoding="utf-8").splitlines() if line.strip()]


def choose_hsp_chain(hsps: list[Hsp]) -> list[Hsp]:
    if not hsps:
        return []
    ordered = sorted(hsps, key=lambda h: (h.qend, h.qstart, -h.bitscore))
    prev: list[int] = [-1] * len(ordered)
    for i in range(len(ordered)):
        for j in range(i - 1, -1, -1):
            if ordered[j].qend < ordered[i].qstart:
                prev[i] = j
                break

    best = [0.0] * len(ordered)
    take = [False] * len(ordered)
    for i, hsp in enumerate(ordered):
        include = hsp.nident + (best[prev[i]] if prev[i] != -1 else 0.0)
        exclude = best[i - 1] if i > 0 else 0.0
        if include > exclude:
            best[i] = include
            take[i] = True
        else:
            best[i] = exclude

    chain: list[Hsp] = []
    i = len(ordered) - 1
    while i >= 0:
        if take[i]:
            chain.append(ordered[i])
            i = prev[i]
        else:
            i -= 1
    chain.reverse()
    return chain


def summarize_chain(chain: list[Hsp], query_aa: str) -> dict[str, Any]:
    if not chain:
        return {
            "refined_match_proportion": 0.0,
            "query_coverage_proportion": 0.0,
            "nident_total": 0,
            "covered_aa": 0,
            "hsp_count": 0,
            "chain_ranges": "",
        }
    qlen = len(query_aa)
    covered = set()
    nident_total = 0
    for h in chain:
        covered.update(range(h.qstart, h.qend + 1))
        nident_total += h.nident
    return {
        "refined_match_proportion": nident_total / qlen if qlen else 0.0,
        "query_coverage_proportion": len(covered) / qlen if qlen else 0.0,
        "nident_total": nident_total,
        "covered_aa": len(covered),
        "hsp_count": len(chain),
        "chain_ranges": "; ".join(
            f"q{h.qstart}-{h.qend}:s{h.sstart}-{h.send}:{h.frames}:nid{h.nident}"
            for h in chain
        ),
    }


def main() -> int:
    args = parse_args()
    output_dir = Path(args.output_dir).expanduser()
    summary = load_json(output_dir / "summary.json")
    gt_refs = build_gt_refs(load_json(Path(args.gt_gene_aa_json).expanduser()))
    subtype_rows = {str(row["accession"]): row for row in load_json(Path(args.subtype_genome_na_json).expanduser())}

    targets = [row for row in summary["records"] if float(row["alignment_score"]) < args.threshold]
    results: list[dict[str, Any]] = []
    per_gene_lines = {
        "NS3": ["NS3 tblastn Refinement", "", f"Threshold: < {args.threshold:.2f}", ""],
        "NS5A_NTD": ["NS5A_NTD tblastn Refinement", "", f"Threshold: < {args.threshold:.2f}", ""],
        "NS5B": ["NS5B tblastn Refinement", "", f"Threshold: < {args.threshold:.2f}", ""],
    }

    for row in targets:
        gene = str(row["gene"])
        genotype = str(row["genotype"])
        accession = str(row["accession"])
        query_aa = gt_refs[(genotype, gene)]
        subject_nt = normalize_nt(str(subtype_rows[accession]["sequence"]))
        tag = f"{gene}_{accession}"
        hsps = run_tblastn(query_aa, subject_nt, tag)
        chain = choose_hsp_chain(hsps)
        summary_row = summarize_chain(chain, query_aa)
        result = {
            "gene": gene,
            "genotype": genotype,
            "subtype": row["subtype"],
            "genotype_name": row["genotype_name"],
            "accession": accession,
            "author_year": row["author_year"],
            "original_match_proportion": float(row["alignment_score"]),
            **summary_row,
            "improvement": summary_row["refined_match_proportion"] - float(row["alignment_score"]),
        }
        results.append(result)

        lines = per_gene_lines[gene]
        lines.append(
            f"> {row['genotype_name']} {accession} original={float(row['alignment_score']):.3f} "
            f"tblastn={summary_row['refined_match_proportion']:.3f} coverage={summary_row['query_coverage_proportion']:.3f}"
        )
        lines.append(summary_row["chain_ranges"] or "no_hsps")
        lines.append("")

    csv_path = output_dir / "tblastn_refined_records_below_0_80.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "gene",
            "genotype",
            "subtype",
            "genotype_name",
            "accession",
            "author_year",
            "original_match_proportion",
            "refined_match_proportion",
            "query_coverage_proportion",
            "nident_total",
            "covered_aa",
            "hsp_count",
            "improvement",
            "chain_ranges",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    comparison_files = {}
    for gene, lines in per_gene_lines.items():
        path = output_dir / f"{gene.lower()}_tblastn_refinement_below_0_80.txt"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        comparison_files[gene] = str(path.resolve())

    payload = {
        "threshold": args.threshold,
        "record_count": len(results),
        "csv": str(csv_path.resolve()),
        "comparison_files": comparison_files,
    }
    (output_dir / "tblastn_refinement_summary.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
