#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$REPO_ROOT/.venv/bin/python}"

usage() {
  cat <<'EOF'
Usage:
  EXCEL_FILE=/path/to/HCV_BlastHits.xlsx FASTA_POOL=/path/to/FASTA scripts/run_ns5a_pipeline.sh

Optional environment variables:
  SHEET_NAME
  OUTPUT_DIR
  REFERENCE_FASTA
  SUBTYPE_JSON
  GT_AA_JSON
  MIN_SEQUENCES
  PYTHON_BIN
EOF
}

EXCEL_FILE="${EXCEL_FILE:-}"
FASTA_POOL="${FASTA_POOL:-}"
SHEET_NAME="${SHEET_NAME:-Ref_summary_20260429 (2)}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO_ROOT/outputs}"
REFERENCE_FASTA="${REFERENCE_FASTA:-$REPO_ROOT/HCV_GT_RefSeqs.fasta}"
SUBTYPE_JSON="${SUBTYPE_JSON:-$REPO_ROOT/HCV_Subtype_Refs_By_Genome_NA.json}"
GT_AA_JSON="${GT_AA_JSON:-$REPO_ROOT/HCV_GT_Refs_By_Gene_AA.json}"
MIN_SEQUENCES="${MIN_SEQUENCES:-10}"
TEMP_ROOT="${TEMP_ROOT:-$REPO_ROOT/temp/$(basename "$0" .sh)}"

if [[ -z "$EXCEL_FILE" || -z "$FASTA_POOL" ]]; then
  usage
  exit 1
fi

MATCHED_TXT="$OUTPUT_DIR/NS5A_matched_fasta_files.txt"
STAGE_DIR="$OUTPUT_DIR/NS5A"
mkdir -p "$TEMP_ROOT"
RUN_TMP="$(mktemp -d "$TEMP_ROOT/run.XXXXXX")"
DISCOVERY_TMP="$RUN_TMP/discovery"
STAGE_ARCHIVE="$(mktemp -d "$TEMP_ROOT/stage_archive.XXXXXX")"

mkdir -p "$OUTPUT_DIR"

cleanup() {
  if [[ -d "$STAGE_DIR" ]]; then
    mv "$STAGE_DIR" "$STAGE_ARCHIVE/NS5A_stage_final" 2>/dev/null || true
  fi
  if [[ -n "${AA_TMP_WORKBOOK:-}" && -f "${AA_TMP_WORKBOOK:-}" ]]; then
    rm -f "$AA_TMP_WORKBOOK"
  fi
}
trap cleanup EXIT

if [[ -d "$STAGE_DIR" ]]; then
  mv "$STAGE_DIR" "$STAGE_ARCHIVE/NS5A_stage_prev"
fi
mkdir -p "$STAGE_DIR"

"$PYTHON_BIN" "$REPO_ROOT/excel-refid-fasta-discovery/scripts/find_refid_fastas.py" \
  --excel-file "$EXCEL_FILE" \
  --fasta-dir "$FASTA_POOL" \
  --output-dir "$DISCOVERY_TMP" \
  --numpatients-column 'Num Pts' \
  --positive-column NS5ACount \
  > "$RUN_TMP/discovery.json"

DISCOVERY_DIR="$(find "$DISCOVERY_TMP" -maxdepth 1 -type d -name 'refid_fasta_*' | head -n 1)"
cp "$DISCOVERY_DIR/matched_fasta_files.txt" "$MATCHED_TXT"

while IFS= read -r src; do
  [[ -n "$src" ]] || continue
  cp "$src" "$STAGE_DIR/"
done < "$MATCHED_TXT"

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns5a_gt_allstudies.py" \
  --excel-file "$EXCEL_FILE" \
  --sheet "$SHEET_NAME" \
  --fasta-dir "$STAGE_DIR" \
  --reference-fasta "$REFERENCE_FASTA" \
  --output-dir "$OUTPUT_DIR" \
  --refid-column RefID \
  --refname-column RefName \
  --numpatients-column 'Num Pts' \
  --positive-column NS5ACount \
  > "$RUN_TMP/gt_allstudies.json"

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns5a_subtype_allstudies_wseqs.py" \
  --combined-workbook "$OUTPUT_DIR/NS5A_GT_AllStudies.xlsx" \
  --fasta-dir "$STAGE_DIR" \
  --subtype-json "$SUBTYPE_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$RUN_TMP/subtype_allstudies.json"

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns5a_subtype_with_gt_aa.py" \
  --subtype-workbook "$OUTPUT_DIR/NS5A_Subtype_AllStudies_WSeqs.xlsx" \
  --fasta-dir "$STAGE_DIR" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$RUN_TMP/subtype_with_gt_aa.json"

AA_TMP_WORKBOOK="$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["output_workbook"])' "$RUN_TMP/subtype_with_gt_aa.json")"

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns5a_completeprofiles_tabspergt.py" \
  --input-workbook "$AA_TMP_WORKBOOK" \
  --output-dir "$OUTPUT_DIR" \
  > "$RUN_TMP/completeprofiles.json"

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns5a_gt_ras_profiles.py" \
  --gt-profile-workbook "$OUTPUT_DIR/NS5A_GT_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$RUN_TMP/gt_ras.json"

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns5a_subtype_ras_profiles_nge10.py" \
  --subtype-profile-workbook "$OUTPUT_DIR/NS5A_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  --min-sequences "$MIN_SEQUENCES" \
  > "$RUN_TMP/subtype_ras.json"

echo "NS5A pipeline complete"
echo "matched_fasta_report=$MATCHED_TXT"
echo "output_dir=$OUTPUT_DIR"
echo "run_tmp=$RUN_TMP"
echo "stage_archive=$STAGE_ARCHIVE"
