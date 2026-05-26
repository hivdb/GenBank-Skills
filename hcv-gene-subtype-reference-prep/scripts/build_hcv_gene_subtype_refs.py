#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile

from Bio import Align


TARGET_GENES = ("NS3", "NS5A_NTD", "NS5B")
NS5A_NTD_AA_LENGTH = 213
NEG_INF = -10**9
MATCH_SCORE = 3
MISMATCH_SCORE = -2
STOP_SCORE = -6
AMBIG_SCORE = -3
FS1_PENALTY = -4
FS2_PENALTY = -5
FS3_PENALTY = -3
CODON_DELETE_PENALTY = -4
GENE_AA_SIGNATURES: dict[tuple[str, str], tuple[str, str]] = {
    ("1", "NS3"): ("API", "VVT"),
    ("1", "NS5A_NTD"): ("SGS", "AEA"),
    ("1", "NS5B"): ("SMS", "PNR"),
    ("2", "NS3"): ("API", "VMT"),
    ("2", "NS5A_NTD"): ("GGS", "AET"),
    ("2", "NS5B"): ("SMS", "PAR"),
    ("3", "NS3"): ("API", "VTT"),
    ("3", "NS5A_NTD"): ("SDD", "AET"),
    ("3", "NS5B"): ("SMS", "PAR"),
    ("4", "NS3"): ("API", "VVT"),
    ("4", "NS5A_NTD"): ("AES", "AES"),
    ("4", "NS5B"): ("SMS", "PAR"),
    ("5", "NS3"): ("API", "VIT"),
    ("5", "NS5A_NTD"): ("DGT", "AET"),
    ("5", "NS5B"): ("SMS", "PAR"),
    ("6", "NS3"): ("API", "VIT"),
    ("6", "NS5A_NTD"): ("ASS", "AET"),
    ("6", "NS5B"): ("SMS", "PAR"),
    ("7", "NS3"): ("API", "VTT"),
    ("7", "NS5A_NTD"): ("AGS", "AET"),
    ("7", "NS5B"): ("SMS", "PAR"),
    ("8", "NS3"): ("SPI", "VVT"),
    ("8", "NS5A_NTD"): ("GNS", "AES"),
    ("8", "NS5B"): ("SMS", "PNR"),
}
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
    exact_matches: int
    informative_sites: int
    ignored_query_positions: int
    aligned_reference_aa: str
    aligned_query_aa: str
    aligned_reference_nt: str
    aligned_query_nt: str
    signature_match_count: int
    signature_mismatch_count: int
    expected_begin_aa: str
    expected_end_aa: str


@dataclass
class FrameshiftRefinement:
    refined_nt: str
    refined_aa: str
    frameshift_events: int
    nt_start: int
    nt_end: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build reusable HCV NS3/NS5A_NTD/NS5B genotype nucleotide/amino-acid "
            "and subtype nucleotide/amino-acid reference FASTA files from FASTA/JSON datasets."
        )
    )
    parser.add_argument("--gt-gene-na-fasta", required=True, help="Path to HCV_GT_RefSeqs.fasta")
    parser.add_argument(
        "--subtype-genome-na-json",
        required=True,
        help="Path to HCV_Subtype_Refs_By_Genome_NA.json",
    )
    parser.add_argument("--output-dir", default="outputs", help="Base output directory")
    parser.add_argument(
        "--enable-frameshift-refinement",
        action="store_true",
        help="Enable the slower frameshift refinement step for sequences with stop codons",
    )
    return parser.parse_args()


def make_job_dir(base_output_dir: Path) -> Path:
    path = base_output_dir / "hcv_gene_subtype_refs"
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


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
    if first_na is None or first_na < 1:
        return None
    return (3 - ((first_na - 1) % 3)) % 3


def translate_nt(nt_sequence: str) -> str:
    aa_chars: list[str] = []
    for idx in range(0, len(nt_sequence), 3):
        codon = nt_sequence[idx : idx + 3]
        if len(codon) < 3:
            break
        aa_chars.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa_chars)


def translate_frame(nt_sequence: str, frame: int) -> tuple[str, str]:
    trimmed = nt_sequence[frame:]
    usable = len(trimmed) - (len(trimmed) % 3)
    coding_nt = trimmed[:usable]
    return coding_nt, translate_nt(coding_nt)


def write_fasta(path: Path, entries: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in entries:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, normalize_nt("".join(chunks))))
            header = line[1:].strip()
            chunks = []
            continue
        chunks.append(line.strip())
    if header is not None:
        records.append((header, normalize_nt("".join(chunks))))
    return records


def parse_gt_fasta_header(header: str) -> tuple[str, str]:
    token = header.split()[0]
    match = re.fullmatch(r"HCV(\d+)([A-Z0-9_]+)", token)
    if not match:
        raise RuntimeError(f"Could not parse genotype FASTA header: {header}")
    return match.group(1), match.group(2)


def build_gt_gene_references(records: list[tuple[str, str]]) -> dict[tuple[str, str], str]:
    refs: dict[tuple[str, str], str] = {}

    for header, sequence in records:
        genotype, gene = parse_gt_fasta_header(header)
        if gene == "NS5A":
            split_index = NS5A_NTD_AA_LENGTH * 3
            if len(sequence) <= split_index:
                raise RuntimeError(
                    f"Genotype {genotype} NS5A reference is too short to split into NTD/CTD: {len(sequence)} nt"
                )
            refs[(genotype, "NS5A_NTD")] = sequence[:split_index]
            refs[(genotype, "NS5A_CTD")] = sequence[split_index:]
        elif gene in {"NS3", "NS5B"}:
            refs[(genotype, gene)] = sequence

    for genotype in map(str, range(1, 9)):
        for gene in TARGET_GENES:
            if (genotype, gene) not in refs:
                raise RuntimeError(f"Missing genotype reference for genotype {genotype} gene {gene}")
    return refs


def excel_cell(value: Any) -> str:
    if value is None:
        return "<c/>"
    if isinstance(value, bool):
        return f'<c t="b"><v>{1 if value else 0}</v></c>'
    if isinstance(value, int | float) and not isinstance(value, bool):
        return f"<c><v>{value}</v></c>"
    return f'<c t="inlineStr"><is><t>{escape(str(value))}</t></is></c>'


def write_alignment_scores_xlsx(path: Path, rows: list[dict[str, Any]]) -> None:
    headers = [
        "gene",
        "genotype",
        "subtype",
        "genotype_name",
        "accession",
        "author_year",
        "alignment_score",
        "frame",
        "query_aa_start",
        "query_aa_end",
        "aa_length",
        "nt_length",
        "exact_matches",
        "informative_sites",
        "ignored_query_positions",
        "signature_match_count",
        "signature_mismatch_count",
        "expected_begin_aa",
        "observed_begin_aa",
        "expected_end_aa",
        "observed_end_aa",
        "has_stop_codon",
        "stop_codon_count",
        "stop_codon_positions",
    ]
    table_rows = [headers]
    table_rows.extend([row.get(header) for header in headers] for row in rows)

    sheet_rows: list[str] = []
    for idx, row_values in enumerate(table_rows, start=1):
        cells = "".join(excel_cell(value) for value in row_values)
        sheet_rows.append(f'<row r="{idx}">{cells}</row>')

    sheet_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f"<sheetData>{''.join(sheet_rows)}</sheetData>"
        "</worksheet>"
    )
    workbook_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets><sheet name="alignment_scores" sheetId="1" r:id="rId1"/></sheets>'
        "</workbook>"
    )
    rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        'Target="worksheets/sheet1.xml"/>'
        "</Relationships>"
    )
    root_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
        'Target="xl/workbook.xml"/>'
        "</Relationships>"
    )
    content_types_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        "</Types>"
    )

    with ZipFile(path, "w", compression=ZIP_DEFLATED) as handle:
        handle.writestr("[Content_Types].xml", content_types_xml)
        handle.writestr("_rels/.rels", root_rels_xml)
        handle.writestr("xl/workbook.xml", workbook_xml)
        handle.writestr("xl/_rels/workbook.xml.rels", rels_xml)
        handle.writestr("xl/worksheets/sheet1.xml", sheet_xml)


def gene_slug(gene: str) -> str:
    return gene.lower()


def stop_codon_positions(sequence: str) -> list[int]:
    return [idx for idx, aa in enumerate(sequence, start=1) if aa == "*"]


def aa_pair_score(ref_aa: str, query_aa: str) -> int:
    if query_aa == "X":
        return AMBIG_SCORE
    if query_aa == "*":
        return STOP_SCORE
    return MATCH_SCORE if ref_aa == query_aa else MISMATCH_SCORE


def aa_signature_stats(
    genotype: str,
    gene: str,
    sequence: str,
) -> tuple[int, int, str, str]:
    expected_begin, expected_end = GENE_AA_SIGNATURES[(genotype, gene)]
    observed_begin = sequence[: len(expected_begin)]
    observed_end = sequence[-len(expected_end) :] if sequence else ""

    matches = 0
    mismatches = 0
    if observed_begin == expected_begin:
        matches += 1
    else:
        mismatches += 1
    if observed_end == expected_end:
        matches += 1
    else:
        mismatches += 1
    return matches, mismatches, expected_begin, expected_end


def is_better_result(candidate: AlignmentResult, incumbent: AlignmentResult | None) -> bool:
    if incumbent is None:
        return True
    if candidate.signature_match_count > incumbent.signature_match_count:
        return True
    if candidate.signature_match_count < incumbent.signature_match_count:
        return False
    if candidate.signature_mismatch_count < incumbent.signature_mismatch_count:
        return True
    if candidate.signature_mismatch_count > incumbent.signature_mismatch_count:
        return False
    candidate_stops = len(stop_codon_positions(candidate.extracted_aa))
    incumbent_stops = len(stop_codon_positions(incumbent.extracted_aa))
    if candidate_stops < incumbent_stops:
        return True
    if candidate_stops > incumbent_stops:
        return False
    if candidate.score > incumbent.score:
        return True
    if candidate.score < incumbent.score:
        return False
    if candidate.informative_sites > incumbent.informative_sites:
        return True
    if candidate.informative_sites < incumbent.informative_sites:
        return False
    if candidate.ignored_query_positions < incumbent.ignored_query_positions:
        return True
    return False


def count_alignment_stats(reference: str, query: str) -> tuple[int, int, int]:
    if len(reference) != len(query):
        raise RuntimeError("Identity scoring requires equal-length sequences")
    matches = 0
    informative = 0
    ignored = 0
    for ref_char, query_char in zip(reference, query):
        if query_char in {"X", "*"}:
            ignored += 1
            continue
        informative += 1
        if ref_char == query_char:
            matches += 1
    return matches, informative, ignored


def sequence_identity(reference: str, query: str) -> tuple[float, int, int, int]:
    matches, informative, ignored = count_alignment_stats(reference, query)
    if informative == 0:
        return 0.0, matches, informative, ignored
    return matches / informative, matches, informative, ignored


def build_aligner(mode: str) -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner(mode=mode)
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    return aligner


def alignment_bounds(blocks: Any) -> tuple[int, int]:
    if len(blocks) == 0:
        raise RuntimeError("Alignment contains no aligned blocks")
    return int(blocks[0][0]), int(blocks[-1][1])


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


def aa_alignment_to_codon_alignment(aligned_reference_aa: str, aligned_query_aa: str, reference_nt: str, query_nt: str) -> tuple[str, str]:
    reference_codons = [reference_nt[idx : idx + 3] for idx in range(0, len(reference_nt), 3)]
    query_codons = [query_nt[idx : idx + 3] for idx in range(0, len(query_nt), 3)]

    ref_index = 0
    query_index = 0
    aligned_reference_nt: list[str] = []
    aligned_query_nt: list[str] = []

    for ref_char, query_char in zip(aligned_reference_aa, aligned_query_aa):
        if ref_char == "-":
            aligned_reference_nt.append("---")
        else:
            aligned_reference_nt.append(reference_codons[ref_index])
            ref_index += 1

        if query_char == "-":
            aligned_query_nt.append("---")
        else:
            aligned_query_nt.append(query_codons[query_index])
            query_index += 1

    return "".join(aligned_reference_nt), "".join(aligned_query_nt)


def discover_gene_window_nt(reference_nt: str, query_nt: str, flank_nt: int = 30) -> tuple[int, int]:
    local_alignment = build_aligner("local").align(reference_nt, query_nt)[0]
    ref_start, ref_end = alignment_bounds(local_alignment.aligned[0])
    query_start, query_end = alignment_bounds(local_alignment.aligned[1])

    window_start = max(0, query_start - ref_start - flank_nt)
    window_end = min(len(query_nt), query_end + (len(reference_nt) - ref_end) + flank_nt)
    return window_start, window_end


def codon_align_gene(
    query_nt: str,
    reference_nt: str,
    reference_aa: str,
) -> tuple[float, str, str, str, str, int, int, int]:
    query_aa = translate_nt(query_nt)
    global_alignment = build_aligner("global").align(reference_aa, query_aa)[0]
    aligned_reference_aa, aligned_query_aa = alignment_strings(reference_aa, query_aa, global_alignment.coordinates)
    aligned_reference_nt, aligned_query_nt = aa_alignment_to_codon_alignment(
        aligned_reference_aa,
        aligned_query_aa,
        reference_nt,
        query_nt,
    )
    score, matches, informative, ignored = sequence_identity(aligned_reference_aa, aligned_query_aa)
    return (
        score,
        aligned_reference_aa,
        aligned_query_aa,
        aligned_reference_nt,
        aligned_query_nt,
        matches,
        informative,
        ignored,
    )


def frameshift_refine(reference_aa: str, nt_window: str) -> FrameshiftRefinement:
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
                aa = CODON_TABLE.get(codon, "X")
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

            if j + 3 <= n:
                score = current + FS3_PENALTY
                if score > dp[i][j + 3]:
                    dp[i][j + 3] = score
                    trace[i][j + 3] = ("skip3", i, j, "")

            score = current + CODON_DELETE_PENALTY
            if score > dp[i + 1][j]:
                dp[i + 1][j] = score
                trace[i + 1][j] = ("delete_ref", i, j, "")

    for j in range(0, n + 1):
        if dp[m][j] > best_score:
            best_score = dp[m][j]
            best_end_j = j

    if best_score <= NEG_INF // 2:
        raise RuntimeError("Frameshift refinement could not align the target window")

    aa_chars: list[str] = []
    kept_nt_chunks: list[str] = []
    frameshift_events = 0
    nt_start = best_end_j
    nt_end = best_end_j
    i = m
    j = best_end_j
    while i > 0 or j > 0:
        step = trace[i][j]
        if step is None:
            break
        op, prev_i, prev_j, aa = step
        if op == "codon":
            aa_chars.append(aa)
            kept_nt_chunks.append(nt_window[prev_j:j])
            nt_start = prev_j
            nt_end = max(nt_end, j)
        else:
            frameshift_events += 1
        i, j = prev_i, prev_j

    aa_chars.reverse()
    kept_nt_chunks.reverse()
    return FrameshiftRefinement(
        refined_nt="".join(kept_nt_chunks),
        refined_aa="".join(aa_chars),
        frameshift_events=frameshift_events,
        nt_start=nt_start,
        nt_end=nt_end,
    )


def build_alignment_result(
    extracted_nt: str,
    reference_nt: str,
    genotype: str,
    gene: str,
    query_aa_start: int,
    frame: int,
) -> AlignmentResult:
    reference_aa = translate_nt(reference_nt)
    extracted_aa = translate_nt(extracted_nt)
    signature_match_count, signature_mismatch_count, expected_begin_aa, expected_end_aa = aa_signature_stats(
        genotype,
        gene,
        extracted_aa,
    )
    (
        score,
        aligned_reference_aa,
        aligned_query_aa,
        aligned_reference_nt,
        aligned_query_nt,
        matches,
        informative,
        ignored,
    ) = codon_align_gene(extracted_nt, reference_nt, reference_aa)
    query_aa_end = query_aa_start + len(extracted_aa)
    return AlignmentResult(
        score=score,
        frame=frame,
        query_aa_start=query_aa_start,
        query_aa_end=query_aa_end,
        extracted_nt=extracted_nt,
        extracted_aa=extracted_aa,
        exact_matches=matches,
        informative_sites=informative,
        ignored_query_positions=ignored,
        aligned_reference_aa=aligned_reference_aa,
        aligned_query_aa=aligned_query_aa,
        aligned_reference_nt=aligned_reference_nt,
        aligned_query_nt=aligned_query_nt,
        signature_match_count=signature_match_count,
        signature_mismatch_count=signature_mismatch_count,
        expected_begin_aa=expected_begin_aa,
        expected_end_aa=expected_end_aa,
    )


def frame_offsets_by_preference(window_start: int, preferred_frame: int | None) -> list[int]:
    preferred_offset = None
    if preferred_frame is not None:
        preferred_offset = (preferred_frame - (window_start % 3)) % 3
    offsets = [preferred_offset] if preferred_offset is not None else []
    offsets.extend(offset for offset in range(3) if offset != preferred_offset)
    return offsets


def extract_gene_from_subtype(
    genotype: str,
    gene: str,
    nt_sequence: str,
    reference_nt: str,
    enable_frameshift_refinement: bool,
    preferred_frame: int | None = None,
) -> AlignmentResult:
    reference_aa = translate_nt(reference_nt)
    window_start, window_end = discover_gene_window_nt(reference_nt, nt_sequence)
    nt_window = nt_sequence[window_start:window_end]
    if not nt_window:
        raise RuntimeError("No nucleotide window could be extracted from subtype sequence")

    candidates: list[AlignmentResult] = []
    last_error: Exception | None = None
    for frame_offset in frame_offsets_by_preference(window_start, preferred_frame):
        coding_nt, _translated = translate_frame(nt_window, frame_offset)
        if not coding_nt:
            continue
        try:
            query_aa_start = (window_start + frame_offset) // 3
            candidates.append(
                build_alignment_result(
                    extracted_nt=coding_nt,
                    reference_nt=reference_nt,
                    genotype=genotype,
                    gene=gene,
                    query_aa_start=query_aa_start,
                    frame=(window_start + frame_offset) % 3,
                )
            )
        except Exception as exc:
            last_error = exc
    if enable_frameshift_refinement:
        try:
            refined = frameshift_refine(reference_aa, nt_window)
            if refined.refined_nt:
                query_aa_start = (window_start + refined.nt_start) // 3
                candidates.append(
                    build_alignment_result(
                        extracted_nt=refined.refined_nt,
                        reference_nt=reference_nt,
                        genotype=genotype,
                        gene=gene,
                        query_aa_start=query_aa_start,
                        frame=(window_start + refined.nt_start) % 3,
                    )
                )
        except Exception as exc:
            last_error = exc

    if not candidates:
        if last_error is not None:
            raise RuntimeError(str(last_error))
        raise RuntimeError("No alignment result could be computed")

    best: AlignmentResult | None = None
    for candidate in candidates:
        if is_better_result(candidate, best):
            best = candidate
    if best is None:
        raise RuntimeError("No alignment result could be selected")
    return best


def build_alignment_marker(reference: str, query: str) -> str:
    chars: list[str] = []
    for ref_char, query_char in zip(reference, query):
        if query_char in {"X", "*"}:
            chars.append(" ")
        elif ref_char == query_char:
            chars.append("|")
        else:
            chars.append(".")
    return "".join(chars)


def write_alignment_report(path: Path, alignments: list[dict[str, Any]], width_codons: int = 68) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for idx, row in enumerate(alignments, start=1):
            header = (
                f"[{idx}] gene={row['gene']} genotype={row['genotype']} subtype={row['subtype']} "
                f"accession={row['accession']} source={row['author_year']}"
            )
            handle.write(header + "\n")
            handle.write(
                "score={score:.6f} frame={frame} aa_range={start}-{end} aa_length={aa_length} "
                "exact_matches={matches}/{informative} ignored_query_positions={ignored} "
                "signature_matches={signature_matches} signature_mismatches={signature_mismatches}\n".format(
                    score=row["alignment_score"],
                    frame=row["frame"],
                    start=row["query_aa_start"],
                    end=row["query_aa_end"],
                    aa_length=row["aa_length"],
                    matches=row["exact_matches"],
                    informative=row["informative_sites"],
                    ignored=row["ignored_query_positions"],
                    signature_matches=row["signature_match_count"],
                    signature_mismatches=row["signature_mismatch_count"],
                )
            )
            handle.write(
                "expected_begin={expected_begin} observed_begin={observed_begin} "
                "expected_end={expected_end} observed_end={observed_end}\n".format(
                    expected_begin=row["expected_begin_aa"],
                    observed_begin=row["observed_begin_aa"],
                    expected_end=row["expected_end_aa"],
                    observed_end=row["observed_end_aa"],
                )
            )
            marker = build_alignment_marker(row["aligned_reference_aa"], row["aligned_query_aa"])
            for start in range(0, len(row["aligned_reference_aa"]), width_codons):
                aa_stop = start + width_codons
                handle.write(f"REF_AA {start + 1:>4} {row['aligned_reference_aa'][start:aa_stop]}\n")
                handle.write(f"MAT_AA {'':>4} {marker[start:aa_stop]}\n")
                handle.write(f"QRY_AA {start + 1:>4} {row['aligned_query_aa'][start:aa_stop]}\n")
                handle.write("\n")
            handle.write("\n")


def main() -> int:
    args = parse_args()
    gt_gene_na_fasta = Path(args.gt_gene_na_fasta).expanduser()
    subtype_genome_na_json = Path(args.subtype_genome_na_json).expanduser()
    base_output_dir = Path(args.output_dir)

    if not gt_gene_na_fasta.exists():
        raise RuntimeError(f"Genotype NA FASTA not found: {gt_gene_na_fasta}")
    if not subtype_genome_na_json.exists():
        raise RuntimeError(f"Subtype genome NA JSON not found: {subtype_genome_na_json}")

    base_output_dir.mkdir(parents=True, exist_ok=True)
    output_dir = make_job_dir(base_output_dir)

    gt_refs = build_gt_gene_references(read_fasta(gt_gene_na_fasta))
    subtype_rows = load_json(subtype_genome_na_json)

    gt_nt_fasta_entries: dict[str, list[tuple[str, str]]] = {gene: [] for gene in TARGET_GENES}
    gt_aa_fasta_entries: dict[str, list[tuple[str, str]]] = {gene: [] for gene in TARGET_GENES}
    for genotype in map(str, range(1, 9)):
        for gene in TARGET_GENES:
            header = f"gene={gene}|genotype={genotype}|source=HCV_GT_RefSeqs.fasta"
            gt_nt_fasta_entries[gene].append((header, gt_refs[(genotype, gene)]))
            gt_aa_fasta_entries[gene].append((header, translate_nt(gt_refs[(genotype, gene)])))
    for gene, entries in gt_nt_fasta_entries.items():
        write_fasta(output_dir / f"hcv_gt_gene_refs_{gene_slug(gene)}_na.fasta", entries)
    for gene, entries in gt_aa_fasta_entries.items():
        write_fasta(output_dir / f"hcv_gt_gene_refs_{gene_slug(gene)}_aa.fasta", entries)

    subtype_nt_entries: dict[str, list[tuple[str, str]]] = {gene: [] for gene in TARGET_GENES}
    subtype_aa_entries: dict[str, list[tuple[str, str]]] = {gene: [] for gene in TARGET_GENES}
    summary_rows: list[dict[str, Any]] = []
    alignment_report_rows: list[dict[str, Any]] = []

    for row in subtype_rows:
        genotype_name = str(row["genotypeName"])
        genotype, subtype = parse_genotype_and_subtype(genotype_name)
        nt_sequence = normalize_nt(str(row["sequence"]))
        accession = str(row.get("accession", ""))
        author_year = str(row.get("authorYear", ""))
        preferred_frame = preferred_frame_from_first_na(row.get("firstNA"))

        for gene in TARGET_GENES:
            reference_nt = gt_refs[(genotype, gene)]
            result = extract_gene_from_subtype(
                genotype=genotype,
                gene=gene,
                nt_sequence=nt_sequence,
                reference_nt=reference_nt,
                enable_frameshift_refinement=args.enable_frameshift_refinement,
                preferred_frame=preferred_frame,
            )
            stop_positions = stop_codon_positions(result.extracted_aa)

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
                    "exact_matches": result.exact_matches,
                    "informative_sites": result.informative_sites,
                    "ignored_query_positions": result.ignored_query_positions,
                    "signature_match_count": result.signature_match_count,
                    "signature_mismatch_count": result.signature_mismatch_count,
                    "expected_begin_aa": result.expected_begin_aa,
                    "observed_begin_aa": result.extracted_aa[: len(result.expected_begin_aa)],
                    "expected_end_aa": result.expected_end_aa,
                    "observed_end_aa": result.extracted_aa[-len(result.expected_end_aa) :],
                    "has_stop_codon": bool(stop_positions),
                    "stop_codon_count": len(stop_positions),
                    "stop_codon_positions": ",".join(map(str, stop_positions)),
                }
            )
            alignment_report_rows.append(
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
                    "exact_matches": result.exact_matches,
                    "informative_sites": result.informative_sites,
                    "ignored_query_positions": result.ignored_query_positions,
                    "signature_match_count": result.signature_match_count,
                    "signature_mismatch_count": result.signature_mismatch_count,
                    "expected_begin_aa": result.expected_begin_aa,
                    "observed_begin_aa": result.extracted_aa[: len(result.expected_begin_aa)],
                    "expected_end_aa": result.expected_end_aa,
                    "observed_end_aa": result.extracted_aa[-len(result.expected_end_aa) :],
                    "has_stop_codon": bool(stop_positions),
                    "stop_codon_count": len(stop_positions),
                    "stop_codon_positions": ",".join(map(str, stop_positions)),
                    "aligned_reference_aa": result.aligned_reference_aa,
                    "aligned_query_aa": result.aligned_query_aa,
                    "aligned_reference_nt": result.aligned_reference_nt,
                    "aligned_query_nt": result.aligned_query_nt,
                }
            )

    for gene, entries in subtype_nt_entries.items():
        write_fasta(output_dir / f"hcv_subtype_gene_refs_{gene_slug(gene)}_na.fasta", entries)
    for gene, entries in subtype_aa_entries.items():
        write_fasta(output_dir / f"hcv_subtype_gene_refs_{gene_slug(gene)}_aa.fasta", entries)
    for gene in TARGET_GENES:
        slug = gene_slug(gene)
        gene_summary_rows = [row for row in summary_rows if row["gene"] == gene]
        gene_alignment_rows = [row for row in alignment_report_rows if row["gene"] == gene]
        write_alignment_scores_xlsx(output_dir / f"{slug}_alignment_scores.xlsx", gene_summary_rows)
        write_alignment_report(output_dir / f"{slug}_alignment_views.txt", gene_alignment_rows)

    summary = {
        "gt_gene_na_fasta": str(gt_gene_na_fasta.resolve()),
        "subtype_genome_na_json": str(subtype_genome_na_json.resolve()),
        "target_genes": list(TARGET_GENES),
        "frameshift_refinement_enabled": bool(args.enable_frameshift_refinement),
        "genotype_reference_count": sum(len(entries) for entries in gt_nt_fasta_entries.values()),
        "subtype_record_count": len(subtype_rows),
        "subtype_gene_sequence_count": sum(len(entries) for entries in subtype_nt_entries.values()),
        "outputs": {
            "gt_gene_refs_na_fastas": {
                gene: str((output_dir / f"hcv_gt_gene_refs_{gene_slug(gene)}_na.fasta").resolve())
                for gene in TARGET_GENES
            },
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
            "alignment_scores_xlsx": {
                gene: str((output_dir / f"{gene_slug(gene)}_alignment_scores.xlsx").resolve())
                for gene in TARGET_GENES
            },
            "alignment_views_txt": {
                gene: str((output_dir / f"{gene_slug(gene)}_alignment_views.txt").resolve())
                for gene in TARGET_GENES
            },
        },
        "records": summary_rows,
    }
    write_json(output_dir / "summary.json", summary)

    print(json.dumps(
        {
            "output_dir": str(output_dir.resolve()),
            "genotype_reference_count": sum(len(entries) for entries in gt_nt_fasta_entries.values()),
            "subtype_gene_sequence_count": sum(len(entries) for entries in subtype_nt_entries.values()),
            "summary_json": str((output_dir / "summary.json").resolve()),
        },
        indent=2,
        ensure_ascii=True,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
