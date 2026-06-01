#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill


RESISTANCE_POSITIONS = [36, 41, 43, 54, 55, 56, 80, 122, 155, 156, 158, 166, 168, 170, 175]
EXCLUDED_AAS = {"X", "*"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build subtype-level NS3 resistance-position AA profile summary in Excel."
    )
    parser.add_argument("--subtype-profile-workbook", required=True)
    parser.add_argument("--gt-aa-json", required=True)
    parser.add_argument("--output-dir", default="outputs")
    return parser.parse_args()


def script_temp_dir() -> Path:
    path = Path("temp") / "hcv-ns3-build-workflow" / Path(__file__).stem
    path.mkdir(parents=True, exist_ok=True)
    return path


def sanitize_label(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in value).strip("._-") or "job"


def make_job_dir(base_output_dir: Path, workbook_path: Path) -> Path:
    label = sanitize_label(f"{workbook_path.stem}_ns3_subtype_resistance_profile")
    job_dir = base_output_dir / label
    if job_dir.exists():
        shutil.rmtree(job_dir)
    job_dir.mkdir(parents=True, exist_ok=True)
    return job_dir


def format_freq(value: float) -> str:
    if value <= 0:
        return "0"
    sig = 2 if value >= 1.0 else 1
    decimals = max(0, sig - 1 - math.floor(math.log10(value)))
    text = f"{round(value, decimals):.{decimals}f}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text


def load_consensus_by_gt(json_path: Path) -> dict[str, str]:
    rows = json.loads(json_path.read_text(encoding="utf-8"))
    consensus: dict[str, str] = {}
    for row in rows:
        name = str(row.get("name", ""))
        match = re.fullmatch(r"HCV([1-8])NS3", name)
        if match:
            consensus[match.group(1)] = str(row.get("refSequence", "")).strip().upper()
    return consensus


def load_subtype_profile_rows(
    workbook_path: Path,
) -> tuple[
    dict[str, dict[str, dict[int, list[tuple[str, float]]]]],
    dict[str, dict[str, int]],
    dict[str, dict[str, dict[int, int]]],
]:
    wb = load_workbook(workbook_path, read_only=True, data_only=True)
    profile_rows: dict[str, dict[str, dict[int, list[tuple[str, float]]]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    subtype_counts: dict[str, dict[str, int]] = defaultdict(dict)
    position_coverage: dict[str, dict[str, dict[int, int]]] = defaultdict(lambda: defaultdict(dict))
    for sheet_name in wb.sheetnames:
        gt = sheet_name.replace("GT", "")
        ws = wb[sheet_name]
        next(ws.iter_rows(values_only=True))
        for row in ws.iter_rows(min_row=2, values_only=True):
            subtype = str(row[0])
            pos = int(row[1])
            denom = int(row[2])
            aa = str(row[3])
            pct = float(row[6])
            if pos in RESISTANCE_POSITIONS:
                profile_rows[gt][subtype][pos].append((aa, pct))
                position_coverage[gt][subtype][pos] = denom
            current = subtype_counts[gt].get(subtype, 0)
            if denom > current:
                subtype_counts[gt][subtype] = denom
    wb.close()
    return profile_rows, subtype_counts, position_coverage


def build_grid(
    consensus_by_gt: dict[str, str],
    profile_rows: dict[str, dict[str, dict[int, list[tuple[str, float]]]]],
    subtype_counts: dict[str, dict[str, int]],
    position_coverage: dict[str, dict[str, dict[int, int]]],
) -> list[list[str]]:
    grid: list[list[str]] = []
    ordered_subtypes: list[tuple[str, str]] = []
    for gt in sorted(subtype_counts, key=int):
        for subtype in sorted(subtype_counts[gt]):
            ordered_subtypes.append((gt, subtype))

    for gt, subtype in ordered_subtypes:
        consensus_seq = consensus_by_gt[gt]
        pos_variants: dict[int, list[str]] = {}
        max_depth = 0
        for pos in RESISTANCE_POSITIONS:
            variants = [
                f"{aa}-{format_freq(pct)}"
                for aa, pct in sorted(profile_rows[gt][subtype].get(pos, []), key=lambda item: (-item[1], item[0]))
                if aa not in EXCLUDED_AAS and pct > 0.1
            ]
            pos_variants[pos] = variants
            max_depth = max(max_depth, len(variants))

        grid.append([f"GT{gt}_{subtype} ({subtype_counts[gt][subtype]})"] + [""] * len(RESISTANCE_POSITIONS))
        grid.append(["Position"] + [str(pos) for pos in RESISTANCE_POSITIONS])
        grid.append(["Reference"] + [consensus_seq[pos - 1] for pos in RESISTANCE_POSITIONS])
        coverage_row = ["Coverage"]
        total_sequences = subtype_counts[gt][subtype]
        for pos in RESISTANCE_POSITIONS:
            covered = position_coverage[gt][subtype].get(pos, 0)
            pct = (100.0 * covered / total_sequences) if total_sequences else 0.0
            coverage_row.append(f"{covered}/{total_sequences} ({pct:.1f}%)")
        grid.append(coverage_row)
        for depth in range(max_depth):
            row = ["Consensus" if depth == 0 else f"Rank{depth + 1}"]
            for pos in RESISTANCE_POSITIONS:
                variants = pos_variants[pos]
                row.append(variants[depth] if depth < len(variants) else "")
            grid.append(row)
        grid.append([""] + [""] * len(RESISTANCE_POSITIONS))
    return grid


def write_excel(path: Path, grid: list[list[str]]) -> None:
    wb = Workbook()
    ws = wb.active
    ws.title = "Subtype_Resistance_Profile"
    block_fill = PatternFill(fill_type="solid", fgColor="D9EAF7")
    header_fill = PatternFill(fill_type="solid", fgColor="F2F2F2")
    consensus_fill = PatternFill(fill_type="solid", fgColor="E2F0D9")
    bold = Font(bold=True)

    for row_idx, row in enumerate(grid, start=1):
        ws.append(row)
        first = row[0]
        if first.startswith("GT"):
            for cell in ws[row_idx]:
                cell.fill = block_fill
                cell.font = bold
        elif first == "Position":
            for cell in ws[row_idx]:
                cell.fill = header_fill
                cell.font = bold
        elif first == "Reference":
            for cell in ws[row_idx]:
                cell.fill = consensus_fill
                cell.font = bold
        for cell in ws[row_idx]:
            cell.alignment = Alignment(horizontal="center")

    from openpyxl.utils import get_column_letter
    for col in range(1, len(RESISTANCE_POSITIONS) + 2):
        ws.column_dimensions[get_column_letter(col)].width = 12 if col > 1 else 18
    wb.save(path)


def main() -> int:
    args = parse_args()
    subtype_profile_workbook = Path(args.subtype_profile_workbook).expanduser()
    gt_aa_json = Path(args.gt_aa_json).expanduser()
    output_dir = Path(args.output_dir)
    script_temp_dir()

    consensus_by_gt = load_consensus_by_gt(gt_aa_json)
    profile_rows, subtype_counts, position_coverage = load_subtype_profile_rows(subtype_profile_workbook)
    grid = build_grid(consensus_by_gt, profile_rows, subtype_counts, position_coverage)

    output_dir.mkdir(parents=True, exist_ok=True)
    excel_path = output_dir / "NS3_Subtype_RAS_Profiles.xlsx"
    write_excel(excel_path, grid)

    summary = {
        "excel": str(excel_path.resolve()),
        "positions": RESISTANCE_POSITIONS,
        "frequency_threshold_percent": 0.1,
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
