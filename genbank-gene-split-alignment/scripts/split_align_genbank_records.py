#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

from Bio import Align
from Bio import Phylo
from Bio import SeqIO
from Bio.Seq import Seq


EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NUCLEOTIDE_PATTERN = re.compile(r"^[ACGTURYKMSWBDHVN]+$", re.IGNORECASE)
AMINO_ACID_PATTERN = re.compile(r"^[ABCDEFGHIKLMNPQRSTVWXYZ*\-]+$", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract nucleotide sequences from GenBank records, align to all reference genes, and split outputs into per-gene FASTA files."
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--accession", action="append", help="GenBank accession; repeat for multiple accessions")
    source.add_argument("--gb-file", help="Local GenBank file path")
    parser.add_argument("--reference-fasta", help="Path to the nucleotide reference FASTA file")
    parser.add_argument("--output-dir", default="outputs", help="Base output directory")
    parser.add_argument("--email", default="", help="Optional contact email for NCBI E-utilities")
    parser.add_argument("--tool", default="genbank-gene-split-alignment", help="Tool name for NCBI E-utilities")
    parser.add_argument("--fully-rerun", action="store_true", help="Ignore cached `.gb` files and recompute everything")
    return parser.parse_args()


def require_reference(args: argparse.Namespace) -> None:
    if not args.reference_fasta:
        raise SystemExit("Missing required input: reference FASTA path")


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


def cohort_name_from_accessions(accessions: list[str]) -> str:
    cleaned = [acc.replace("/", "_") for acc in accessions]
    job_id = make_job_id("|".join(cleaned))
    if len(cleaned) <= 3:
        label = "_".join(cleaned)
    else:
        label = f"{cleaned[0]}_to_{cleaned[-1]}_{len(cleaned)}records"
    return f"gene_split_{sanitize_label(label)}_{job_id}"


def cohort_name_from_gb_file(gb_file: Path) -> str:
    label = sanitize_label(gb_file.stem)
    job_id = make_job_id(str(gb_file.resolve()))
    return f"gene_split_{label}_{job_id}"


def load_json(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def find_existing_cohort_dir_for_gb_file(base_output_dir: Path, gb_file: Path) -> Path | None:
    target = str(gb_file.resolve())
    for summary_path in sorted(base_output_dir.glob("cohort_*/summary.json"), reverse=True):
        payload = load_json(summary_path)
        if not payload:
            continue
        if payload.get("source_gb_file") == target:
            return summary_path.parent
    return None


def find_existing_cohort_dir_for_accessions(base_output_dir: Path, accessions: list[str]) -> Path | None:
    target = list(accessions)
    for summary_path in sorted(base_output_dir.glob("cohort_*/summary.json"), reverse=True):
        payload = load_json(summary_path)
        if not payload:
            continue
        if payload.get("source_accessions") == target:
            return summary_path.parent
    return None


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
            pass
        else:
            gb_path.write_text(fetch_genbank_record(accession, email, tool), encoding="utf-8")
        records.extend(list(SeqIO.parse(gb_path, "genbank")))
    return records


def load_records_from_gb_file(gb_file: Path) -> list[Any]:
    if gb_file.is_dir():
        records: list[Any] = []
        for path in sorted(gb_file.glob("*.gb*")):
            records.extend(list(SeqIO.parse(path, "genbank")))
        return records
    return list(SeqIO.parse(gb_file, "genbank"))


def read_reference_genes(reference_fasta: Path) -> list[dict[str, str]]:
    genes: list[dict[str, str]] = []
    for record in SeqIO.parse(reference_fasta, "fasta"):
        sequence = str(record.seq).upper().replace("-", "")
        if not sequence:
            raise RuntimeError(f"Reference FASTA entry '{record.id}' is empty")
        if NUCLEOTIDE_PATTERN.match(sequence):
            sequence_type = "nucleotide"
        elif AMINO_ACID_PATTERN.match(sequence):
            sequence_type = "amino_acid"
        else:
            raise RuntimeError(
                f"Reference FASTA entry '{record.id}' is not recognized as nucleotide or amino-acid sequence"
            )
        genes.append(
            {
                "gene": record.id,
                "description": record.description,
                "sequence": sequence,
                "sequence_type": sequence_type,
            }
        )
    if not genes:
        raise RuntimeError(f"No reference gene sequences were found in FASTA: {reference_fasta}")
    return genes


def make_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    return aligner


def compute_nucleotide_alignment(reference_seq: str, query_seq: str) -> tuple[dict[str, Any], str]:
    aligner = make_aligner()
    alignment = aligner.align(reference_seq, query_seq)[0]
    ref_segments = alignment.aligned[0]
    query_segments = alignment.aligned[1]

    ref_start = int(ref_segments[0][0]) + 1
    ref_end = int(ref_segments[-1][1])
    query_start = int(query_segments[0][0]) + 1
    query_end = int(query_segments[-1][1])

    projected = ["-"] * len(reference_seq)
    for (r0, r1), (q0, q1) in zip(ref_segments, query_segments):
        projected[r0:r1] = list(query_seq[q0:q1])

    summary = {
        "score": alignment.score,
        "reference_start": ref_start,
        "reference_end": ref_end,
        "query_start": query_start,
        "query_end": query_end,
        "reference_aligned_length": sum(end - start for start, end in ref_segments),
        "query_aligned_length": sum(end - start for start, end in query_segments),
    }
    return summary, "".join(projected)


def translate_query_frames(query_seq: str) -> list[dict[str, Any]]:
    frames: list[dict[str, Any]] = []
    clean_seq = query_seq.upper().replace("-", "")
    for frame in range(3):
        coding_len = ((len(clean_seq) - frame) // 3) * 3
        if coding_len <= 0:
            continue
        frame_nt = clean_seq[frame : frame + coding_len]
        aa_seq = str(Seq(frame_nt).translate(to_stop=False))
        frames.append(
            {
                "frame": frame,
                "nt_sequence": frame_nt,
                "aa_sequence": aa_seq,
            }
        )
    return frames


def project_codons_from_aa_alignment(reference_aa: str, query_nt: str) -> tuple[dict[str, Any], str]:
    best_summary = None
    best_projected = ""
    best_query_nt = ""
    aligner = make_aligner()
    translated_frames = translate_query_frames(query_nt)
    for frame_info in translated_frames:
        aa_alignment = aligner.align(reference_aa, frame_info["aa_sequence"])[0]
        ref_segments = aa_alignment.aligned[0]
        query_segments = aa_alignment.aligned[1]
        if len(ref_segments) == 0 or len(query_segments) == 0:
            continue

        projected = ["---"] * len(reference_aa)
        for (r0, r1), (q0, q1) in zip(ref_segments, query_segments):
            for ref_pos, query_pos in zip(range(r0, r1), range(q0, q1)):
                codon = frame_info["nt_sequence"][query_pos * 3 : query_pos * 3 + 3]
                if len(codon) == 3:
                    projected[ref_pos] = codon

        summary = {
            "score": aa_alignment.score,
            "reference_start": int(ref_segments[0][0]) * 3 + 1,
            "reference_end": int(ref_segments[-1][1]) * 3,
            "query_start": frame_info["frame"] + int(query_segments[0][0]) * 3 + 1,
            "query_end": frame_info["frame"] + int(query_segments[-1][1]) * 3,
            "reference_aligned_length": sum(end - start for start, end in ref_segments) * 3,
            "query_aligned_length": sum(end - start for start, end in query_segments) * 3,
            "query_frame": frame_info["frame"],
            "alignment_mode": "codon_projected_from_aa_reference",
        }
        if best_summary is None or summary["score"] > best_summary["score"]:
            best_summary = summary
            best_projected = "".join(projected)
            best_query_nt = frame_info["nt_sequence"]

    if best_summary is None:
        raise RuntimeError("No codon alignment could be computed against amino-acid reference")
    return best_summary, best_projected


def best_gene_alignment(record, reference_genes: list[dict[str, str]]) -> tuple[dict[str, Any], str]:
    query_seq = str(record.seq).upper().replace("-", "")
    best_summary = None
    best_projected = ""
    best_gene = None
    for gene in reference_genes:
        if gene["sequence_type"] == "nucleotide":
            summary, projected = compute_nucleotide_alignment(gene["sequence"], query_seq)
            summary["alignment_mode"] = "nucleotide_to_nucleotide"
        else:
            summary, projected = project_codons_from_aa_alignment(gene["sequence"], query_seq)
        if best_summary is None or summary["score"] > best_summary["score"]:
            best_summary = summary
            best_projected = projected
            best_gene = gene
    if best_summary is None or best_gene is None:
        raise RuntimeError(f"No alignment result could be computed for {record.id}")
    result = {
        "accession": accession_from_record(record),
        "version": record.id,
        "matched_gene": best_gene["gene"],
        "reference_description": best_gene["description"],
        "reference_sequence_type": best_gene["sequence_type"],
        **best_summary,
    }
    return result, best_projected


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def write_gene_fastas(path: Path, grouped: dict[str, list[tuple[str, str]]]) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for gene, entries in grouped.items():
        gene_path = path / f"{gene}.fasta"
        with gene_path.open("w", encoding="utf-8") as fh:
            for accession, aligned_sequence in entries:
                fh.write(f">{accession}\n{aligned_sequence}\n")


def find_fasttree_binary() -> str | None:
    for candidate in ("FastTree", "fasttree", "FastTreeMP"):
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    return None


def run_fasttree(gene_fasta: Path, tree_path: Path) -> None:
    binary = find_fasttree_binary()
    if not binary:
        raise SystemExit("FastTree was not found on PATH. Please install FastTree and rerun this workflow.")
    with tree_path.open("w", encoding="utf-8") as tree_fh:
        subprocess.run(
            [binary, "-nt", "-gtr", "-boot", "1000", str(gene_fasta)],
            check=True,
            stdout=tree_fh,
            stderr=subprocess.PIPE,
            text=True,
        )


def render_tree_figure(tree_path: Path, png_path: Path, svg_path: Path | None = None) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tree = Phylo.read(tree_path, "newick")
    tip_count = len(tree.get_terminals())
    height = max(4.0, min(0.18 * tip_count + 2.0, 40.0))
    fig = plt.figure(figsize=(12.0, height))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_axis_off()
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=True)
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    if svg_path is not None:
        fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)


def build_trees_for_genes(alignments_dir: Path, trees_dir: Path, grouped: dict[str, list[tuple[str, str]]]) -> dict[str, dict[str, str]]:
    trees_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, dict[str, str]] = {}
    for gene, entries in sorted(grouped.items()):
        if len(entries) < 2:
            outputs[gene] = {
                "status": "skipped",
                "reason": "need at least two sequences",
            }
            continue
        gene_fasta = alignments_dir / f"{gene}.fasta"
        tree_path = trees_dir / f"{gene}.treefile"
        png_path = trees_dir / f"{gene}.png"
        svg_path = trees_dir / f"{gene}.svg"
        run_fasttree(gene_fasta, tree_path)
        render_tree_figure(tree_path, png_path, svg_path)
        outputs[gene] = {
            "status": "created",
            "treefile": str(tree_path),
            "figure_png": str(png_path),
            "figure_svg": str(svg_path),
        }
    return outputs


def main() -> int:
    args = parse_args()
    require_reference(args)

    base_output_dir = Path(args.output_dir)
    base_output_dir.mkdir(parents=True, exist_ok=True)
    reference_genes = read_reference_genes(Path(args.reference_fasta).expanduser())

    if args.accession:
        cohort_dir = find_existing_cohort_dir_for_accessions(base_output_dir, args.accession)
        if cohort_dir is None:
            cohort_name = cohort_name_from_accessions(args.accession)
            cohort_dir = base_output_dir / cohort_name
            cohort_dir.mkdir(parents=True, exist_ok=True)
        else:
            cohort_name = cohort_dir.name
        records = load_records_from_accessions(args.accession, cohort_dir, args.email, args.tool, args.fully_rerun)
    else:
        gb_file = Path(args.gb_file).expanduser()
        cohort_dir = find_existing_cohort_dir_for_gb_file(base_output_dir, gb_file)
        if cohort_dir is None:
            cohort_name = cohort_name_from_gb_file(gb_file)
            cohort_dir = base_output_dir / cohort_name
            cohort_dir.mkdir(parents=True, exist_ok=True)
        else:
            cohort_name = cohort_dir.name
        records = load_records_from_gb_file(gb_file)

    if not records:
        raise RuntimeError("No GenBank records were found to process")

    rows: list[dict[str, Any]] = []
    grouped: dict[str, list[tuple[str, str]]] = defaultdict(list)
    gene_counts: Counter[str] = Counter()
    for record in records:
        result, projected = best_gene_alignment(record, reference_genes)
        rows.append(result)
        grouped[result["matched_gene"]].append((result["accession"], projected))
        gene_counts[result["matched_gene"]] += 1

    summary_csv = cohort_dir / "gene_alignment_summary.csv"
    write_csv(
        summary_csv,
        rows,
        [
            "accession",
            "version",
            "matched_gene",
            "reference_description",
            "reference_sequence_type",
            "reference_start",
            "reference_end",
            "query_start",
            "query_end",
            "reference_aligned_length",
            "query_aligned_length",
            "score",
            "query_frame",
            "alignment_mode",
        ],
    )

    alignments_dir = cohort_dir / "alignments"
    write_gene_fastas(alignments_dir, grouped)
    trees_dir = cohort_dir / "trees"
    tree_outputs = build_trees_for_genes(alignments_dir, trees_dir, grouped)

    summary = {
        "cohort_name": cohort_name,
        "record_count": len(rows),
        "gene_counts": dict(gene_counts),
        "gene_alignment_summary_csv": str(summary_csv),
        "alignments_dir": str(alignments_dir),
        "trees_dir": str(trees_dir),
        "tree_outputs": tree_outputs,
        "input_mode": "accessions" if args.accession else "gb_file",
    }
    write_json(cohort_dir / "gene_alignment_summary.json", summary)

    print(json.dumps(summary, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
