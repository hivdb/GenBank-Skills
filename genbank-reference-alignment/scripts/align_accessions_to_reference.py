#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

from Bio import Align
from Bio import SeqIO


EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Align one or more GenBank accessions to all sequences in a reference FASTA and pick the best match."
    )
    query_group = parser.add_mutually_exclusive_group(required=True)
    query_group.add_argument("--accession", action="append", help="GenBank accession; repeat for multiple accessions")
    query_group.add_argument("--query-fasta", help="FASTA file containing one or more query sequences to align")
    parser.add_argument("--reference-fasta", help="Path to the reference FASTA file")
    parser.add_argument("--gene", help="Optional gene filter to limit which reference FASTA headers are considered")
    parser.add_argument("--output-dir", default="outputs", help="Base output directory")
    parser.add_argument("--email", default="", help="Optional contact email to send to NCBI E-utilities")
    parser.add_argument("--tool", default="genbank-reference-alignment", help="Tool name to send to NCBI E-utilities")
    parser.add_argument("--fully-rerun", action="store_true", help="Ignore cached accession files and recompute everything")
    return parser.parse_args()


def require_user_inputs(args: argparse.Namespace) -> None:
    missing: list[str] = []
    if not args.reference_fasta:
        missing.append("reference FASTA path")
    if missing:
        joined = " and ".join(missing)
        raise SystemExit(f"Missing required input: {joined}")


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


def fetch_fasta(accession: str, email: str, tool: str) -> str:
    query = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "tool": tool,
    }
    if email:
        query["email"] = email
    text = fetch_text(f"{EFETCH_URL}?{urlencode(query)}")
    if not text.strip():
        raise RuntimeError(f"Empty FASTA response for {accession}")
    return text


def read_query_fasta(query_fasta: Path) -> list[dict[str, str]]:
    queries: list[dict[str, str]] = []
    for record in SeqIO.parse(query_fasta, "fasta"):
        queries.append(
            {
                "query_id": record.id,
                "header": record.description,
                "sequence": str(record.seq).upper(),
                "input_type": "fasta",
            }
        )
    if not queries:
        raise RuntimeError(f"No query sequences were found in FASTA: {query_fasta}")
    return queries


def read_reference_candidates(reference_fasta: Path, gene: str | None) -> list[dict[str, str]]:
    candidates: list[dict[str, str]] = []
    gene_lower = gene.lower() if gene else ""
    for record in SeqIO.parse(reference_fasta, "fasta"):
        header = f"{record.id} {record.description}".strip()
        if gene and gene_lower not in header.lower():
            continue
        candidates.append(
            {
                "id": record.id,
                "description": record.description,
                "matched_gene_filter": gene or "",
                "sequence": str(record.seq).upper(),
            }
        )
    if not candidates:
        if gene:
            raise RuntimeError(f"Gene '{gene}' was not found in reference FASTA: {reference_fasta}")
        raise RuntimeError(f"No reference sequences were found in FASTA: {reference_fasta}")
    return candidates


def parse_fasta_sequence(text: str) -> dict[str, str]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise RuntimeError("Invalid FASTA payload")
    header = lines[0][1:]
    sequence = "".join(lines[1:]).upper()
    accession = header.split()[0]
    return {"header": header, "accession": accession, "sequence": sequence}


def compute_alignment(reference_seq: str, query_seq: str) -> dict[str, Any]:
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(reference_seq, query_seq)[0]
    ref_segments = alignment.aligned[0]
    query_segments = alignment.aligned[1]

    ref_start = int(ref_segments[0][0]) + 1
    ref_end = int(ref_segments[-1][1])
    query_start = int(query_segments[0][0]) + 1
    query_end = int(query_segments[-1][1])

    return {
        "score": alignment.score,
        "reference_start": ref_start,
        "reference_end": ref_end,
        "query_start": query_start,
        "query_end": query_end,
        "reference_aligned_length": sum(end - start for start, end in ref_segments),
        "query_aligned_length": sum(end - start for start, end in query_segments),
    }


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=True, default=json_default) + "\n",
        encoding="utf-8",
    )


def json_default(value: Any) -> Any:
    if hasattr(value, "item"):
        return value.item()
    raise TypeError(f"Object of type {value.__class__.__name__} is not JSON serializable")


def build_single_report(result: dict[str, Any]) -> str:
    return "\n".join(
        [
            "GenBank Reference Alignment",
            "",
            f"Query ID: {result['query_id']}",
            f"Input type: {result['input_type']}",
            f"Accession: {result.get('accession')}",
            f"Matched reference gene: {result['matched_reference_gene']}",
            f"Reference record: {result['reference_record_id']}",
            f"Reference description: {result['reference_description']}",
            f"Reference range: {result['reference_start']}-{result['reference_end']}",
            f"Query range: {result['query_start']}-{result['query_end']}",
            f"Alignment score: {result['score']}",
            f"Reference aligned length: {result['reference_aligned_length']}",
            f"Query aligned length: {result['query_aligned_length']}",
        ]
    ) + "\n"


def load_cached_alignment(output_dir: Path) -> dict[str, Any] | None:
    alignment_path = output_dir / "alignment.json"
    if not alignment_path.exists():
        return None
    cached = json.loads(alignment_path.read_text(encoding="utf-8"))
    required_keys = {
        "query_id",
        "matched_reference_gene",
        "reference_record_id",
        "reference_start",
        "reference_end",
        "query_start",
        "query_end",
        "score",
    }
    if not required_keys.issubset(cached):
        return None
    return cached


def process_query(
    query: dict[str, str],
    base_output_dir: Path,
    reference_candidates: list[dict[str, str]],
    email: str,
    tool: str,
    fully_rerun: bool,
) -> dict[str, Any]:
    query_id = query["query_id"]
    query_dir = base_output_dir / query_id.replace("/", "_")
    query_dir.mkdir(parents=True, exist_ok=True)

    cached = None if fully_rerun else load_cached_alignment(query_dir)
    if cached is not None:
        return cached

    fasta_path = query_dir / "sequence.fasta"
    fasta: dict[str, str]

    if query["input_type"] == "accession":
        accession = query["accession"]
        record_path = query_dir / "record.gb"
        if record_path.exists() and not fully_rerun:
            record_text = record_path.read_text(encoding="utf-8")
        else:
            record_text = fetch_genbank_record(accession, email, tool)
            record_path.write_text(record_text, encoding="utf-8")

        if fasta_path.exists() and not fully_rerun:
            fasta_text = fasta_path.read_text(encoding="utf-8")
        else:
            fasta_text = fetch_fasta(accession, email, tool)
            fasta_path.write_text(fasta_text, encoding="utf-8")
        fasta = parse_fasta_sequence(fasta_text)
    else:
        fasta = {
            "header": query["header"],
            "accession": query_id,
            "sequence": query["sequence"],
        }
        if fully_rerun or not fasta_path.exists():
            fasta_path.write_text(f">{query['header']}\n{query['sequence']}\n", encoding="utf-8")

    best_reference = None
    best_alignment = None
    for reference_candidate in reference_candidates:
        alignment = compute_alignment(reference_candidate["sequence"], fasta["sequence"])
        if best_alignment is None or alignment["score"] > best_alignment["score"]:
            best_reference = reference_candidate
            best_alignment = alignment

    if best_reference is None or best_alignment is None:
        raise RuntimeError(f"No alignment result could be computed for {query_id}")

    result = {
        "query_id": query_id,
        "input_type": query["input_type"],
        "accession": query.get("accession"),
        "fasta_header": fasta["header"],
        "matched_reference_gene": best_reference["id"],
        "reference_record_id": best_reference["id"],
        "reference_description": best_reference["description"],
        "reference_gene_filter": best_reference["matched_gene_filter"],
        **best_alignment,
    }
    write_json(query_dir / "alignment.json", result)
    return result


def write_multi_csv(path: Path, results: list[dict[str, Any]]) -> None:
    fieldnames = [
        "query_id",
        "input_type",
        "accession",
        "matched_reference_gene",
        "reference_record_id",
        "reference_description",
        "reference_start",
        "reference_end",
        "query_start",
        "query_end",
        "reference_aligned_length",
        "query_aligned_length",
        "score",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow({key: result.get(key) for key in fieldnames})


def build_batch_csv_name(identifiers: list[str], source_label: str = "") -> str:
    cleaned = [identifier.replace("/", "_") for identifier in identifiers]
    if len(cleaned) <= 3:
        suffix = "_".join(cleaned)
    else:
        suffix = f"{cleaned[0]}_to_{cleaned[-1]}_{len(cleaned)}accessions"
    if source_label:
        return f"alignment_batch_{source_label}_{suffix}.csv"
    return f"alignment_batch_{suffix}.csv"


def main() -> int:
    args = parse_args()
    require_user_inputs(args)

    base_output_dir = Path(args.output_dir)
    base_output_dir.mkdir(parents=True, exist_ok=True)

    reference_candidates = read_reference_candidates(Path(args.reference_fasta), args.gene)
    queries: list[dict[str, str]]
    source_label = ""
    if args.query_fasta:
        query_fasta = Path(args.query_fasta)
        queries = read_query_fasta(query_fasta)
        source_label = query_fasta.stem
    else:
        queries = [
            {
                "query_id": accession,
                "header": accession,
                "sequence": "",
                "input_type": "accession",
                "accession": accession,
            }
            for accession in args.accession or []
        ]

    results = [
        process_query(
            query,
            base_output_dir,
            reference_candidates,
            args.email,
            args.tool,
            args.fully_rerun,
        )
        for query in queries
    ]

    if len(results) == 1:
        report_path = base_output_dir / results[0]["query_id"].replace("/", "_") / "alignment_report.txt"
        report_path.write_text(build_single_report(results[0]), encoding="utf-8")
        print(
            json.dumps(
                {
                    "mode": "single",
                    "report": str(report_path),
                    "alignment": str(base_output_dir / results[0]["query_id"].replace("/", "_") / "alignment.json"),
                },
                indent=2,
                ensure_ascii=True,
            )
        )
        return 0

    csv_path = base_output_dir / build_batch_csv_name([result["query_id"] for result in results], source_label)
    write_multi_csv(csv_path, results)
    print(
        json.dumps(
            {
                "mode": "multi",
                "csv": str(csv_path),
                "message": f"Look for the batch CSV file at {csv_path}",
                "count": len(results),
            },
            indent=2,
            ensure_ascii=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
