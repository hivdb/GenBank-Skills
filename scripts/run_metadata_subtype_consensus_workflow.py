#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from Bio import Align
from openpyxl import Workbook


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA_CSV = REPO_ROOT / "temp/build_accessions_metadata_csv/Accessions_metadata.csv"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "temp/metadata_subtype_consensus_workflow"
DEFAULT_SUBTYPE_JSON = REPO_ROOT / "HCV_Subtype_Refs_By_Genome_NA.json"
DEFAULT_GT_AA_JSON = REPO_ROOT / "HCV_GT_Refs_By_Gene_AA.json"
DEFAULT_GT_AA_FASTA = REPO_ROOT / "HCV_GT_Refs_NS3_NS5A_NTD_NS5B_AA.fasta"
DEFAULT_REFERENCE_FASTA = REPO_ROOT / "HCV_GT_RefSeqs.fasta"
GENES = ("NS3", "NS5A", "NS5B")
GT_REF_GENE_BY_GENE = {"NS3": "NS3", "NS5A": "NS5A_NTD", "NS5B": "NS5B"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build subtype complete-profile workbooks and subtype consensus FASTA files "
            "from an accession metadata CSV and a RefID-organized FASTA directory."
        )
    )
    parser.add_argument("--fasta-dir", help="Directory containing FASTA files named with RefID prefixes")
    parser.add_argument("--metadata-csv", default=str(DEFAULT_METADATA_CSV))
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--subtype-json", default=str(DEFAULT_SUBTYPE_JSON))
    parser.add_argument("--gt-aa-json", default=str(DEFAULT_GT_AA_JSON))
    parser.add_argument("--gt-aa-fasta", default=str(DEFAULT_GT_AA_FASTA))
    parser.add_argument("--reference-fasta", default=str(DEFAULT_REFERENCE_FASTA))
    parser.add_argument("--python-bin", default=sys.executable)
    parser.add_argument("--genes", nargs="+", choices=GENES, default=list(GENES))
    parser.add_argument("--min-aligned-nt", type=int, default=200)
    parser.add_argument("--min-aa-overlap", type=int, default=80)
    parser.add_argument("--alignment-width", type=int, default=80)
    parser.add_argument(
        "--only-consensus-alignments",
        action="store_true",
        help="Only align existing *_Subtype_Consensus.fasta files to genotype AA references and write text reports.",
    )
    return parser.parse_args()


def is_gene_present(value: str | None) -> bool:
    return (value or "").strip().lower() in {"yes", "y", "true", "1", "present"}


def load_metadata_accessions(metadata_csv: Path, gene: str) -> tuple[dict[str, set[str]], dict[str, int]]:
    accessions_by_refid: dict[str, set[str]] = {}
    seen: set[tuple[str, str]] = set()
    skipped_no_gene = 0
    skipped_missing_id = 0
    skipped_duplicate = 0

    with metadata_csv.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"RefID", "Accession", gene}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise RuntimeError(f"{metadata_csv} is missing columns: {', '.join(missing)}")

        for row in reader:
            if not is_gene_present(row.get(gene)):
                skipped_no_gene += 1
                continue
            refid = (row.get("RefID") or "").strip()
            accession = (row.get("Accession") or "").strip()
            if not refid or not accession:
                skipped_missing_id += 1
                continue
            key = (refid, accession)
            if key in seen:
                skipped_duplicate += 1
                continue
            seen.add(key)
            accessions_by_refid.setdefault(refid, set()).add(accession)

    return accessions_by_refid, {
        "selected_accessions": sum(len(accessions) for accessions in accessions_by_refid.values()),
        "selected_refids": len(accessions_by_refid),
        "skipped_no_gene": skipped_no_gene,
        "skipped_missing_refid_or_accession": skipped_missing_id,
        "skipped_duplicate": skipped_duplicate,
    }


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks).upper()))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append("".join(line.split()))
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    return records


def accession_from_header(header: str) -> str:
    return header.split()[0]


def write_fasta(path: Path, entries: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in entries:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")


def refid_sort_key(refid: str) -> tuple[int, str]:
    return (int(refid) if refid.isdigit() else 10**12, refid)


def build_refid_to_fasta(fasta_dir: Path) -> dict[str, Path]:
    mapping: dict[str, Path] = {}
    for path in sorted(fasta_dir.rglob("*")):
        if not path.is_file() or "_" not in path.name:
            continue
        refid = path.name.split("_", 1)[0]
        mapping.setdefault(refid, path)
    return mapping


def stage_gene_fastas(
    fasta_dir: Path,
    accessions_by_refid: dict[str, set[str]],
    stage_dir: Path,
    gene: str,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    if stage_dir.exists():
        shutil.rmtree(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)
    refid_to_fasta = build_refid_to_fasta(fasta_dir)
    seed_rows: list[dict[str, Any]] = []
    missing_fasta = 0
    missing_sequence = 0

    for refid in sorted(accessions_by_refid, key=refid_sort_key):
        fasta_path = refid_to_fasta.get(refid)
        wanted = accessions_by_refid[refid]
        if fasta_path is None:
            missing_fasta += len(wanted)
            continue
        entries = [(header, seq) for header, seq in parse_fasta(fasta_path) if accession_from_header(header) in wanted]
        found = {accession_from_header(header) for header, _seq in entries}
        missing_sequence += len(wanted - found)
        if not entries:
            continue
        staged_path = stage_dir / f"{refid}_metadata_{gene}.fasta"
        write_fasta(staged_path, entries)
        seed_rows.append(
            {
                "RefID": refid,
                "RefName": f"metadata_{refid}",
                "NumPts": len(entries),
                f"{gene}Count": len(entries),
            }
        )
    return seed_rows, {
        "staged_refids": len(seed_rows),
        "staged_accessions": sum(int(row["NumPts"]) for row in seed_rows),
        "missing_fasta_accessions": missing_fasta,
        "missing_sequence_accessions": missing_sequence,
    }


def write_seed_workbook(path: Path, rows: list[dict[str, Any]], gene: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    wb = Workbook()
    ws = wb.active
    ws.title = "Metadata_GT"
    header = ["RefID", "RefName", "NumPts", f"{gene}Count"]
    ws.append(header)
    for row in rows:
        ws.append([row.get(field, "") for field in header])
    wb.save(path)


def run_json_command(command: list[str], summary_path: Path) -> dict[str, Any]:
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    summary_path.write_text(completed.stdout.strip() + "\n", encoding="utf-8")
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"Command did not emit JSON: {' '.join(command)}") from exc


def workflow_scripts(gene: str) -> dict[str, Path]:
    lower = gene.lower()
    return {
        "gt": REPO_ROOT / f"scripts/build_{lower}_gt_allstudies.py",
        "subtype": REPO_ROOT / f"scripts/build_{lower}_subtype_allstudies_wseqs.py",
        "aa": REPO_ROOT / f"scripts/build_{lower}_subtype_with_gt_aa.py",
        "profiles": REPO_ROOT / f"scripts/build_{lower}_completeprofiles_tabspergt.py",
    }


def build_aa_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    return aligner


def alignment_strings(reference: str, query: str, coordinates: Any) -> tuple[str, str]:
    ref_chunks: list[str] = []
    query_chunks: list[str] = []
    ref_points = coordinates[0]
    query_points = coordinates[1]
    for idx in range(len(ref_points) - 1):
        ref_start = int(ref_points[idx])
        ref_end = int(ref_points[idx + 1])
        query_start = int(query_points[idx])
        query_end = int(query_points[idx + 1])
        ref_span = ref_end - ref_start
        query_span = query_end - query_start
        if ref_span > 0 and query_span > 0:
            if ref_span != query_span:
                raise RuntimeError("Unexpected alignment segment with unequal spans")
            ref_chunks.append(reference[ref_start:ref_end])
            query_chunks.append(query[query_start:query_end])
        elif ref_span > 0:
            ref_chunks.append(reference[ref_start:ref_end])
            query_chunks.append("-" * ref_span)
        elif query_span > 0:
            ref_chunks.append("-" * query_span)
            query_chunks.append(query[query_start:query_end])
    return "".join(ref_chunks), "".join(query_chunks)


def alignment_marker(reference: str, query: str) -> str:
    chars: list[str] = []
    for ref_char, query_char in zip(reference, query):
        if ref_char == "-" or query_char == "-" or query_char in {"X", "*"}:
            chars.append(" ")
        elif ref_char == query_char:
            chars.append("|")
        else:
            chars.append(".")
    return "".join(chars)


def parse_consensus_genotype(header: str) -> str:
    match = re.match(r"GT([1-8])_", header)
    if not match:
        raise RuntimeError(f"Could not parse genotype from subtype consensus header: {header}")
    return match.group(1)


def parse_gt_ref_header(header: str) -> tuple[str, str] | None:
    parts = [part.strip() for part in header.split("|")]
    if len(parts) < 2:
        return None
    match = re.fullmatch(r"genotype\s+([1-8])", parts[0], flags=re.IGNORECASE)
    if not match:
        return None
    return match.group(1), parts[1]


def load_gt_aa_fasta_refs(path: Path) -> dict[tuple[str, str], str]:
    refs: dict[tuple[str, str], str] = {}
    for header, sequence in parse_fasta(path):
        parsed = parse_gt_ref_header(header)
        if parsed is None:
            continue
        refs[parsed] = sequence
    return refs


def write_consensus_alignment_report(
    gene: str,
    consensus_fasta: Path,
    gt_refs: dict[tuple[str, str], str],
    gt_aa_fasta: Path,
    output_path: Path,
    width: int,
) -> dict[str, Any]:
    ref_gene = GT_REF_GENE_BY_GENE[gene]
    records = parse_fasta(consensus_fasta)
    aligner = build_aa_aligner()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    aligned_count = 0

    with output_path.open("w", encoding="utf-8") as handle:
        for idx, (header, sequence) in enumerate(records, start=1):
            gt = parse_consensus_genotype(header)
            reference = gt_refs.get((gt, ref_gene))
            if reference is None:
                raise RuntimeError(f"Missing genotype {gt} {ref_gene} AA reference in {gt_aa_fasta}")

            alignment = aligner.align(reference, sequence)[0]
            aligned_ref, aligned_query = alignment_strings(reference, sequence, alignment.coordinates)
            marker = alignment_marker(aligned_ref, aligned_query)
            informative = sum(
                1
                for ref_char, query_char in zip(aligned_ref, aligned_query)
                if ref_char != "-" and query_char not in {"-", "X", "*"}
            )
            matches = marker.count("|")
            identity = matches / informative if informative else 0.0

            handle.write(f"[{idx}] gene={gene} subtype_consensus={header} genotype=GT{gt} reference={ref_gene}\n")
            handle.write(
                f"score={alignment.score:.3f} identity={matches}/{informative} ({identity:.3%}) "
                f"reference_length={len(reference)} consensus_length={len(sequence)}\n"
            )
            for start in range(0, len(aligned_ref), width):
                end = start + width
                handle.write(f"REF  {start + 1:>5} {aligned_ref[start:end]}\n")
                handle.write(f"           {marker[start:end]}\n")
                handle.write(f"CONS {start + 1:>5} {aligned_query[start:end]}\n\n")
            aligned_count += 1

    return {
        "gene": gene,
        "consensus_fasta": str(consensus_fasta.resolve()),
        "alignment_txt": str(output_path.resolve()),
        "records_aligned": aligned_count,
    }


def run_consensus_alignment_report(args: argparse.Namespace, output_root: Path, gene: str) -> dict[str, Any]:
    gt_aa_fasta = Path(args.gt_aa_fasta).expanduser()
    gt_refs = load_gt_aa_fasta_refs(gt_aa_fasta)
    gene_dir = output_root / gene
    consensus_fasta = gene_dir / f"{gene}_Subtype_Consensus.fasta"
    if not consensus_fasta.exists():
        raise RuntimeError(f"Missing subtype consensus FASTA: {consensus_fasta}")
    output_path = gene_dir / f"{gene}_Subtype_Consensus_GT_Alignment.txt"
    return write_consensus_alignment_report(
        gene=gene,
        consensus_fasta=consensus_fasta,
        gt_refs=gt_refs,
        gt_aa_fasta=gt_aa_fasta,
        output_path=output_path,
        width=args.alignment_width,
    )


def run_consensus_alignment_reports(args: argparse.Namespace, output_root: Path) -> list[dict[str, Any]]:
    return [run_consensus_alignment_report(args, output_root, gene) for gene in args.genes]


def run_gene(args: argparse.Namespace, gene: str, output_root: Path) -> dict[str, Any]:
    fasta_dir = Path(args.fasta_dir).expanduser()
    metadata_csv = Path(args.metadata_csv).expanduser()
    subtype_json = Path(args.subtype_json).expanduser()
    gt_aa_json = Path(args.gt_aa_json).expanduser()
    reference_fasta = Path(args.reference_fasta).expanduser()
    gene_dir = output_root / gene
    gene_dir.mkdir(parents=True, exist_ok=True)

    accessions_by_refid, metadata_summary = load_metadata_accessions(metadata_csv, gene)
    stage_dir = gene_dir / "fasta_stage"
    seed_rows, stage_summary = stage_gene_fastas(fasta_dir, accessions_by_refid, stage_dir, gene)
    seed_workbook = gene_dir / f"{gene}_Metadata_GT_Seed.xlsx"
    write_seed_workbook(seed_workbook, seed_rows, gene)
    if not seed_rows:
        return {
            "gene": gene,
            "seed_workbook": str(seed_workbook.resolve()),
            "metadata": metadata_summary,
            "fasta_stage": stage_summary,
            "status": "skipped_no_rows",
        }

    scripts = workflow_scripts(gene)
    gt_summary = run_json_command(
        [
            args.python_bin,
            str(scripts["gt"]),
            "--excel-file",
            str(seed_workbook),
            "--sheet",
            "Metadata_GT",
            "--fasta-dir",
            str(stage_dir),
            "--reference-fasta",
            str(reference_fasta),
            "--output-dir",
            str(gene_dir),
            "--refid-column",
            "RefID",
            "--refname-column",
            "RefName",
            "--numpatients-column",
            "NumPts",
            "--min-aligned-nt",
            str(args.min_aligned_nt),
        ],
        gene_dir / f"{gene}_gt_assignment_summary.json",
    )

    subtype_summary = run_json_command(
        [
            args.python_bin,
            str(scripts["subtype"]),
            "--combined-workbook",
            str(gt_summary["combined_xlsx"]),
            "--fasta-dir",
            str(stage_dir),
            "--subtype-json",
            str(subtype_json),
            "--output-dir",
            str(gene_dir),
            "--min-aligned-nt",
            str(args.min_aligned_nt),
        ],
        gene_dir / f"{gene}_subtype_assignment_summary.json",
    )

    aa_workbook = gene_dir / f"{gene}_Subtype_With_GT_AA.xlsx"
    aa_summary = run_json_command(
        [
            args.python_bin,
            str(scripts["aa"]),
            "--subtype-workbook",
            str(subtype_summary["output_workbook"]),
            "--fasta-dir",
            str(fasta_dir),
            "--gt-aa-json",
            str(gt_aa_json),
            "--output-dir",
            str(gene_dir),
            "--output-workbook",
            str(aa_workbook),
            "--min-aa-overlap",
            str(args.min_aa_overlap),
        ],
        gene_dir / f"{gene}_aa_extraction_summary.json",
    )

    profile_summary = run_json_command(
        [
            args.python_bin,
            str(scripts["profiles"]),
            "--input-workbook",
            str(aa_workbook),
            "--output-dir",
            str(gene_dir),
        ],
        gene_dir / f"{gene}_completeprofiles_summary.json",
    )

    consensus_summary = run_json_command(
        [
            args.python_bin,
            str(REPO_ROOT / "scripts/export_subtype_consensus_fasta.py"),
            "--gene",
            gene,
            "--subtype-profile-workbook",
            str(profile_summary["subtype_workbook"]),
            "--output-dir",
            str(gene_dir),
        ],
        gene_dir / f"{gene}_subtype_consensus_summary.json",
    )
    consensus_alignment_summary = run_consensus_alignment_report(args, output_root, gene)

    return {
        "gene": gene,
        "status": "complete",
        "seed_workbook": str(seed_workbook.resolve()),
        "metadata": metadata_summary,
        "fasta_stage": stage_summary,
        "gt_assignment": gt_summary,
        "subtype_assignment": subtype_summary,
        "aa_extraction": aa_summary,
        "complete_profiles": profile_summary,
        "subtype_consensus": consensus_summary,
        "subtype_consensus_alignment": consensus_alignment_summary,
    }


def main() -> int:
    args = parse_args()
    if not args.only_consensus_alignments and not args.fasta_dir:
        raise SystemExit("--fasta-dir is required unless --only-consensus-alignments is used")
    output_root = Path(args.output_root).expanduser()
    output_root.mkdir(parents=True, exist_ok=True)

    if args.only_consensus_alignments:
        alignment_summaries = run_consensus_alignment_reports(args, output_root)
        summary = {
            "output_root": str(output_root.resolve()),
            "gt_aa_fasta": str(Path(args.gt_aa_fasta).expanduser().resolve()),
            "genes": alignment_summaries,
        }
        summary_path = output_root / "consensus_alignment_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(summary, indent=2))
        return 0

    summaries = [run_gene(args, gene, output_root) for gene in args.genes]
    summary = {
        "metadata_csv": str(Path(args.metadata_csv).expanduser().resolve()),
        "fasta_dir": str(Path(args.fasta_dir or "").expanduser().resolve()),
        "output_root": str(output_root.resolve()),
        "genes": summaries,
    }
    summary_path = output_root / "workflow_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
