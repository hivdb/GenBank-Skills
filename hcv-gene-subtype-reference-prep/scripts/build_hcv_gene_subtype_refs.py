#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


TARGET_GENES = ("NS3", "NS5A_NTD", "NS5B")
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


@dataclass
class AlignmentResult:
    score: float
    frame: int
    query_aa_start: int
    query_aa_end: int
    extracted_nt: str
    extracted_aa: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build reusable HCV NS3/NS5A_NTD/NS5B genotype amino-acid and subtype "
            "nucleotide/amino-acid reference FASTA files from JSON datasets."
        )
    )
    parser.add_argument("--gt-gene-aa-json", required=True, help="Path to HCV_GT_Refs_By_Gene_AA.json")
    parser.add_argument(
        "--subtype-genome-na-json",
        required=True,
        help="Path to HCV_Subtype_Refs_By_Genome_NA.json",
    )
    parser.add_argument("--output-dir", default="outputs", help="Base output directory")
    return parser.parse_args()


def sanitize_label(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return text.strip("._-") or "job"


def make_job_dir(base_output_dir: Path) -> Path:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    path = base_output_dir / f"hcv_gene_subtype_refs_{timestamp}"
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_gt_from_name(name: str) -> str:
    match = re.match(r"HCV(\d+)", name)
    if not match:
        raise RuntimeError(f"Could not parse genotype from name: {name}")
    return match.group(1)


def parse_genotype_and_subtype(genotype_name: str) -> tuple[str, str]:
    match = re.match(r"Genotype(\d+)([A-Za-z0-9]+)$", genotype_name)
    if not match:
        raise RuntimeError(f"Could not parse genotype/subtype from genotypeName: {genotype_name}")
    genotype = match.group(1)
    subtype = f"{genotype}{match.group(2)}"
    return genotype, subtype


def normalize_nt(sequence: str) -> str:
    return re.sub(r"[^ACGTUNacgtun]", "", sequence).upper().replace("U", "T")


def preferred_frame_from_first_na(first_na: int | None) -> int | None:
    if first_na is None:
        return None
    if first_na < 1:
        return None
    return (3 - ((first_na - 1) % 3)) % 3


def translate_frame(nt_sequence: str, frame: int) -> tuple[str, str]:
    trimmed = nt_sequence[frame:]
    usable = len(trimmed) - (len(trimmed) % 3)
    coding_nt = trimmed[:usable]
    aa_chars: list[str] = []
    for idx in range(0, usable, 3):
        aa_chars.append(CODON_TABLE.get(coding_nt[idx : idx + 3], "X"))
    return coding_nt, "".join(aa_chars)


def write_fasta(path: Path, entries: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in entries:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def gene_slug(gene: str) -> str:
    return gene.lower()


def sequence_identity(reference: str, query: str) -> float:
    if len(reference) != len(query):
        raise RuntimeError("Identity scoring requires equal-length sequences")
    if not reference:
        return 0.0
    matches = 0
    informative = 0
    for ref_char, query_char in zip(reference, query):
        if query_char in {"X", "*"}:
            continue
        informative += 1
        if ref_char == query_char:
            matches += 1
    if informative == 0:
        return 0.0
    return matches / informative


def choose_seed_length(reference_length: int) -> int:
    if reference_length >= 600:
        return 12
    if reference_length >= 300:
        return 10
    return 8


def seed_lengths(reference_length: int) -> list[int]:
    primary = choose_seed_length(reference_length)
    candidates = [primary, primary - 2, primary - 4, 6, 5, 4]
    result: list[int] = []
    for value in candidates:
        if value >= 4 and value not in result:
            result.append(value)
    return result


def find_seed_hits(reference: str, query: str, seed_length: int) -> list[tuple[int, int]]:
    hits: list[tuple[int, int]] = []
    for ref_pos in range(0, len(reference) - seed_length + 1):
        seed = reference[ref_pos : ref_pos + seed_length]
        query_pos = query.find(seed)
        while query_pos != -1:
            hits.append((ref_pos, query_pos))
            query_pos = query.find(seed, query_pos + 1)
    return hits


def anchored_window_search(reference: str, query: str) -> tuple[float, int, int]:
    best_score = -1.0
    best_start = -1

    for seed_length in seed_lengths(len(reference)):
        hits = find_seed_hits(reference, query, seed_length)
        if not hits:
            continue
        for ref_pos, query_pos in hits:
            start = query_pos - ref_pos
            end = start + len(reference)
            if start < 0 or end > len(query):
                continue
            window = query[start:end]
            score = sequence_identity(reference, window)
            if score > best_score:
                best_score = score
                best_start = start
        if best_start >= 0:
            return best_score, best_start, best_start + len(reference)

    raise RuntimeError("No exact amino-acid seed hit found between reference and query")


def build_gt_gene_references(rows: list[dict[str, Any]]) -> dict[tuple[str, str], str]:
    refs: dict[tuple[str, str], str] = {}

    for row in rows:
        genotype = parse_gt_from_name(str(row["name"]))
        gene = str(row["abstractGene"])
        sequence = str(row["refSequence"]).strip().upper()

        if gene in {"NS3", "NS5A_NTD", "NS5B"}:
            refs[(genotype, gene)] = sequence

    for genotype in map(str, range(1, 9)):
        for gene in TARGET_GENES:
            if (genotype, gene) not in refs:
                raise RuntimeError(f"Missing genotype reference for genotype {genotype} gene {gene}")
    return refs


def extract_gene_from_subtype(
    nt_sequence: str,
    reference_aa: str,
    preferred_frame: int | None = None,
) -> AlignmentResult:
    best: AlignmentResult | None = None
    frames = [preferred_frame] if preferred_frame is not None else []
    frames.extend(frame for frame in range(3) if frame != preferred_frame)
    last_error = None
    for frame in frames:
        coding_nt, translated = translate_frame(nt_sequence, frame)
        try:
            score, aa_start, aa_end = anchored_window_search(reference_aa, translated)
        except RuntimeError as exc:
            last_error = exc
            continue
        nt_start = aa_start * 3
        nt_end = aa_end * 3
        extracted_nt = coding_nt[nt_start:nt_end]
        extracted_aa = translated[aa_start:aa_end]
        candidate = AlignmentResult(
            score=score,
            frame=frame,
            query_aa_start=aa_start,
            query_aa_end=aa_end,
            extracted_nt=extracted_nt,
            extracted_aa=extracted_aa,
        )
        if best is None or candidate.score > best.score:
            best = candidate

    if best is None:
        if last_error is not None:
            raise RuntimeError(str(last_error))
        raise RuntimeError("No alignment result could be computed")
    return best


def main() -> int:
    args = parse_args()
    gt_gene_aa_json = Path(args.gt_gene_aa_json).expanduser()
    subtype_genome_na_json = Path(args.subtype_genome_na_json).expanduser()
    base_output_dir = Path(args.output_dir)

    if not gt_gene_aa_json.exists():
        raise RuntimeError(f"Genotype AA JSON not found: {gt_gene_aa_json}")
    if not subtype_genome_na_json.exists():
        raise RuntimeError(f"Subtype genome NA JSON not found: {subtype_genome_na_json}")

    base_output_dir.mkdir(parents=True, exist_ok=True)
    output_dir = make_job_dir(base_output_dir)

    gt_rows = load_json(gt_gene_aa_json)
    subtype_rows = load_json(subtype_genome_na_json)

    gt_refs = build_gt_gene_references(gt_rows)

    gt_fasta_entries: dict[str, list[tuple[str, str]]] = {gene: [] for gene in TARGET_GENES}
    for genotype in map(str, range(1, 9)):
        for gene in TARGET_GENES:
            header = f"gene={gene}|genotype={genotype}|source=HCV_GT_Refs_By_Gene_AA"
            gt_fasta_entries[gene].append((header, gt_refs[(genotype, gene)]))
    for gene, entries in gt_fasta_entries.items():
        write_fasta(output_dir / f"hcv_gt_gene_refs_{gene_slug(gene)}_aa.fasta", entries)

    subtype_nt_entries: dict[str, list[tuple[str, str]]] = {gene: [] for gene in TARGET_GENES}
    subtype_aa_entries: dict[str, list[tuple[str, str]]] = {gene: [] for gene in TARGET_GENES}
    summary_rows: list[dict[str, Any]] = []

    for row in subtype_rows:
        genotype_name = str(row["genotypeName"])
        genotype, subtype = parse_genotype_and_subtype(genotype_name)
        nt_sequence = normalize_nt(str(row["sequence"]))
        accession = str(row.get("accession", ""))
        author_year = str(row.get("authorYear", ""))
        preferred_frame = preferred_frame_from_first_na(row.get("firstNA"))

        for gene in TARGET_GENES:
            reference_aa = gt_refs[(genotype, gene)]
            result = extract_gene_from_subtype(nt_sequence, reference_aa, preferred_frame=preferred_frame)

            base_header = (
                f"gene={gene}|genotype={genotype}|subtype={subtype}|"
                f"accession={accession}|genotypeName={genotype_name}|source={author_year}"
            )
            subtype_nt_entries[gene].append((base_header, result.extracted_nt))
            subtype_aa_entries[gene].append((base_header, result.extracted_aa))
            summary_rows.append(
                {
                    "gene": gene,
                    "genotype": genotype,
                    "subtype": subtype,
                    "genotype_name": genotype_name,
                    "accession": accession,
                    "author_year": author_year,
                    "alignment_score": result.score,
                    "frame": result.frame,
                    "query_aa_start": result.query_aa_start + 1,
                    "query_aa_end": result.query_aa_end,
                    "nt_length": len(result.extracted_nt),
                    "aa_length": len(result.extracted_aa),
                }
            )

    for gene, entries in subtype_nt_entries.items():
        write_fasta(output_dir / f"hcv_subtype_gene_refs_{gene_slug(gene)}_na.fasta", entries)
    for gene, entries in subtype_aa_entries.items():
        write_fasta(output_dir / f"hcv_subtype_gene_refs_{gene_slug(gene)}_aa.fasta", entries)

    summary = {
        "gt_gene_aa_json": str(gt_gene_aa_json.resolve()),
        "subtype_genome_na_json": str(subtype_genome_na_json.resolve()),
        "target_genes": list(TARGET_GENES),
        "genotype_reference_count": sum(len(entries) for entries in gt_fasta_entries.values()),
        "subtype_record_count": len(subtype_rows),
        "subtype_gene_sequence_count": sum(len(entries) for entries in subtype_nt_entries.values()),
        "outputs": {
            "gt_gene_refs_aa_fastas": {
                gene: str((output_dir / f"hcv_gt_gene_refs_{gene_slug(gene)}_aa.fasta").resolve())
                for gene in TARGET_GENES
            },
            "subtype_gene_refs_na_fastas": {
                gene: str((output_dir / f"hcv_subtype_gene_refs_{gene_slug(gene)}_na.fasta").resolve())
                for gene in TARGET_GENES
            },
            "subtype_gene_refs_aa_fastas": {
                gene: str((output_dir / f"hcv_subtype_gene_refs_{gene_slug(gene)}_aa.fasta").resolve())
                for gene in TARGET_GENES
            },
        },
        "records": summary_rows,
    }
    write_json(output_dir / "summary.json", summary)

    print(json.dumps(
        {
            "output_dir": str(output_dir.resolve()),
            "genotype_reference_count": sum(len(entries) for entries in gt_fasta_entries.values()),
            "subtype_gene_sequence_count": sum(len(entries) for entries in subtype_nt_entries.values()),
            "summary_json": str((output_dir / "summary.json").resolve()),
        },
        indent=2,
        ensure_ascii=True,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
