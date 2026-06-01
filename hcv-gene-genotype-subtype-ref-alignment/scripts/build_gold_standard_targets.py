#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any


OUTFMT = "6 qseqid sseqid qlen qstart qend sstart send length nident mismatch gaps pident evalue bitscore qseq sseq frames"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build final gold-standard HCV target AA FASTAs from original and tblastn-refined candidates."
    )
    parser.add_argument("--output-dir", required=True, help="Output directory from build_hcv_gene_subtype_refs.py")
    parser.add_argument("--gt-gene-aa-json", required=True, help="Path to HCV_GT_Refs_By_Gene_AA.json")
    parser.add_argument("--subtype-genome-na-json", required=True, help="Path to HCV_Subtype_Refs_By_Genome_NA.json")
    parser.add_argument("--original-accept-threshold", type=float, default=0.80)
    parser.add_argument("--refined-accept-threshold", type=float, default=0.80)
    parser.add_argument("--rescued-accept-threshold", type=float, default=0.70)
    parser.add_argument("--min-improvement", type=float, default=0.15)
    return parser.parse_args()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def normalize_nt(sequence: str) -> str:
    return "".join(base for base in sequence.upper().replace("U", "T") if base in {"A", "C", "G", "T", "N"})


def parse_gt_from_name(name: str) -> str:
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


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks)
            header = line[1:]
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        records[header] = "".join(chunks)
    return records


def parse_header(header: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for part in header.split("|"):
        if "=" in part:
            key, value = part.split("=", 1)
            result[key] = value
    return result


def write_fasta(path: Path, entries: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in entries:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start:start + 70] + "\n")


def run_tblastn(query_aa: str, subject_nt: str, tag: str) -> list[dict[str, Any]]:
    with tempfile.TemporaryDirectory(prefix=f"tblastn_{tag}_") as tmpdir:
        tmp = Path(tmpdir)
        query_path = tmp / "query.faa"
        subject_path = tmp / "subject.fna"
        out_path = tmp / "hits.tsv"
        write_fasta(query_path, [(tag, query_aa)])
        write_fasta(subject_path, [(tag, subject_nt)])
        subprocess.run(
            [
                "tblastn", "-query", str(query_path), "-subject", str(subject_path),
                "-seg", "no", "-comp_based_stats", "F", "-max_hsps", "50", "-evalue", "1e-3",
                "-outfmt", OUTFMT, "-out", str(out_path),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            text=True,
        )
        hits = []
        for line in out_path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            p = line.split("\t")
            hits.append({
                "qlen": int(p[2]), "qstart": int(p[3]), "qend": int(p[4]),
                "sstart": int(p[5]), "send": int(p[6]), "length": int(p[7]),
                "nident": int(p[8]), "bitscore": float(p[13]), "qseq": p[14], "sseq": p[15], "frames": p[16],
            })
        return hits


def choose_hsp_chain(hsps: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not hsps:
        return []
    ordered = sorted(hsps, key=lambda h: (h["qend"], h["qstart"], -h["bitscore"]))
    prev = [-1] * len(ordered)
    for i in range(len(ordered)):
        for j in range(i - 1, -1, -1):
            if ordered[j]["qend"] < ordered[i]["qstart"]:
                prev[i] = j
                break
    best = [0.0] * len(ordered)
    take = [False] * len(ordered)
    for i, h in enumerate(ordered):
        include = h["nident"] + (best[prev[i]] if prev[i] != -1 else 0.0)
        exclude = best[i - 1] if i > 0 else 0.0
        if include > exclude:
            best[i] = include
            take[i] = True
        else:
            best[i] = exclude
    chain = []
    i = len(ordered) - 1
    while i >= 0:
        if take[i]:
            chain.append(ordered[i])
            i = prev[i]
        else:
            i -= 1
    return list(reversed(chain))


def reconstruct_target(query_aa: str, chain: list[dict[str, Any]]) -> tuple[str, float, float, str]:
    qlen = len(query_aa)
    target = ["X"] * qlen
    match_count = 0
    covered = set()
    for h in chain:
        qseq = h["qseq"]
        sseq = h["sseq"]
        qpos = h["qstart"] - 1
        for qa, sa in zip(qseq, sseq):
            if qa != "-":
                if sa != "-":
                    target[qpos] = sa
                    covered.add(qpos)
                    if qa == sa:
                        match_count += 1
                qpos += 1
    target_seq = "".join(target)
    return target_seq, match_count / qlen if qlen else 0.0, len(covered) / qlen if qlen else 0.0, "; ".join(
        f"q{h['qstart']}-{h['qend']}:s{h['sstart']}-{h['send']}:{h['frames']}:nid{h['nident']}" for h in chain
    )


def main() -> int:
    args = parse_args()
    outdir = Path(args.output_dir).expanduser()
    summary = load_json(outdir / "summary.json")
    gt_refs = build_gt_refs(load_json(Path(args.gt_gene_aa_json).expanduser()))
    subtype_rows = {str(r["accession"]): r for r in load_json(Path(args.subtype_genome_na_json).expanduser())}

    original_aa = {}
    for name in ["ns3", "ns5a_ntd", "ns5b"]:
        for header, seq in read_fasta(outdir / f"hcv_subtype_gene_refs_{name}_aa.fasta").items():
            meta = parse_header(header)
            original_aa[(meta["gene"], meta["accession"])] = (header, seq)

    decisions = []
    accepted_by_gene = {"NS3": [], "NS5A_NTD": [], "NS5B": []}
    review = []

    for row in summary["records"]:
        gene = str(row["gene"])
        genotype = str(row["genotype"])
        accession = str(row["accession"])
        original_score = float(row["alignment_score"])
        header, original_seq = original_aa[(gene, accession)]

        method = "original"
        final_seq = original_seq
        final_score = original_score
        coverage = 1.0
        chain_ranges = ""
        status = "accepted_original" if original_score >= args.original_accept_threshold else "manual_review"

        if original_score < args.original_accept_threshold:
            query_aa = gt_refs[(genotype, gene)]
            subject_nt = normalize_nt(str(subtype_rows[accession]["sequence"]))
            chain = choose_hsp_chain(run_tblastn(query_aa, subject_nt, f"{gene}_{accession}"))
            refined_seq, refined_score, coverage, chain_ranges = reconstruct_target(query_aa, chain)
            if refined_score >= args.refined_accept_threshold or (
                refined_score >= args.rescued_accept_threshold and refined_score - original_score >= args.min_improvement
            ):
                method = "tblastn_refined"
                final_seq = refined_seq
                final_score = refined_score
                status = "accepted_tblastn_refined"
            else:
                final_score = refined_score

        record = {
            "gene": gene,
            "genotype": genotype,
            "subtype": row["subtype"],
            "genotype_name": row["genotype_name"],
            "accession": accession,
            "author_year": row["author_year"],
            "status": status,
            "method": method,
            "original_match_proportion": original_score,
            "final_match_proportion": final_score,
            "query_coverage_proportion": coverage,
            "chain_ranges": chain_ranges,
        }
        decisions.append(record)

        if status.startswith("accepted"):
            accepted_by_gene[gene].append((header + f"|selection={method}|final_match={final_score:.3f}", final_seq))
        else:
            review.append(record)

    for gene, slug in [("NS3", "ns3"), ("NS5A_NTD", "ns5a_ntd"), ("NS5B", "ns5b")]:
        accepted_by_gene[gene].sort(key=lambda item: item[0])
        write_fasta(outdir / f"gold_standard_{slug}_aa.fasta", accepted_by_gene[gene])

    fieldnames = [
        "gene", "genotype", "subtype", "genotype_name", "accession", "author_year",
        "status", "method", "original_match_proportion", "final_match_proportion",
        "query_coverage_proportion", "chain_ranges",
    ]
    with (outdir / "gold_standard_selection_table.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in decisions:
            writer.writerow(row)

    with (outdir / "gold_standard_manual_review.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in review:
            writer.writerow(row)

    payload = {
        "accepted_original": sum(1 for r in decisions if r["status"] == "accepted_original"),
        "accepted_tblastn_refined": sum(1 for r in decisions if r["status"] == "accepted_tblastn_refined"),
        "manual_review": len(review),
        "selection_table": str((outdir / "gold_standard_selection_table.csv").resolve()),
        "manual_review_table": str((outdir / "gold_standard_manual_review.csv").resolve()),
        "gold_standard_fastas": {
            "NS3": str((outdir / "gold_standard_ns3_aa.fasta").resolve()),
            "NS5A_NTD": str((outdir / "gold_standard_ns5a_ntd_aa.fasta").resolve()),
            "NS5B": str((outdir / "gold_standard_ns5b_aa.fasta").resolve()),
        },
    }
    (outdir / "gold_standard_summary.json").write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
