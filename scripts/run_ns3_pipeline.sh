#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$REPO_ROOT/.venv/bin/python}"
CONFIG_PATH="$REPO_ROOT/pipeline.local.toml"

if [[ -f "$CONFIG_PATH" ]]; then
  eval "$(python3 "$REPO_ROOT/scripts/load_pipeline_defaults.py" ns3 "$REPO_ROOT")"
fi

usage() {
  cat <<'EOF'
Usage:
  EXCEL_FILE=/path/to/HCV_BlastHits.xlsx FASTA_POOL=/path/to/FASTA [GENBANK_DIR=/path/to/genbank_seq_files] scripts/run_ns3_pipeline.sh

Optional environment variables:
  SHEET_NAME
  OUTPUT_DIR
  REFERENCE_FASTA
  SUBTYPE_JSON
  GT_AA_JSON
  PYTHON_BIN
  TEMP_ROOT

Defaults can also be provided in pipeline.local.toml.
EOF
}

EXCEL_FILE="${EXCEL_FILE:-}"
FASTA_POOL="${FASTA_POOL:-}"
GENBANK_DIR="${GENBANK_DIR:-}"
SHEET_NAME="${SHEET_NAME:-Ref_summary_20260429 (2)}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO_ROOT/outputs}"
REFERENCE_FASTA="${REFERENCE_FASTA:-$REPO_ROOT/HCV_GT_RefSeqs.fasta}"
SUBTYPE_JSON="${SUBTYPE_JSON:-$REPO_ROOT/HCV_Subtype_Refs_By_Genome_NA.json}"
GT_AA_JSON="${GT_AA_JSON:-$REPO_ROOT/HCV_GT_Refs_By_Gene_AA.json}"
TEMP_ROOT="${TEMP_ROOT:-$REPO_ROOT/temp/$(basename "$0" .sh)}"

if [[ -z "$EXCEL_FILE" || -z "$FASTA_POOL" ]]; then
  usage
  exit 1
fi

MATCHED_TXT="$OUTPUT_DIR/NS3_matched_fasta_files.txt"
STAGE_DIR="$TEMP_ROOT/NS3_stage"
DISCOVERY_TMP="$REPO_ROOT/temp/find_refid_fastas"
DISCOVERY_JSON="$DISCOVERY_TMP/discovery_ns3.json"
GT_ALLSTUDIES_JSON="$REPO_ROOT/temp/build_ns3_gt_allstudies/last_run_summary.json"
SOURCEFEATURES_JSON="$REPO_ROOT/temp/build_ns3_sourcefeatures_csv/last_run_summary.json"
SOURCEFEATURES_GROUPED_JSON="$REPO_ROOT/temp/build_ns3_sourcefeatures_grouped_csv/last_run_summary.json"
SUBTYPE_ALLSTUDIES_JSON="$REPO_ROOT/temp/build_ns3_subtype_allstudies_wseqs/last_run_summary.json"
SUBTYPE_WITH_GT_AA_JSON="$REPO_ROOT/temp/build_ns3_subtype_with_gt_aa/last_run_summary.json"
COMPLETEPROFILES_JSON="$REPO_ROOT/temp/build_ns3_completeprofiles_tabspergt/last_run_summary.json"
GT_RAS_JSON="$REPO_ROOT/temp/build_ns3_gt_ras_profiles/last_run_summary.json"
SUBTYPE_RAS_JSON="$REPO_ROOT/temp/build_ns3_subtype_ras_profiles/last_run_summary.json"
mkdir -p "$TEMP_ROOT"
mkdir -p "$DISCOVERY_TMP"
mkdir -p "$(dirname "$GT_ALLSTUDIES_JSON")" "$(dirname "$SOURCEFEATURES_JSON")" "$(dirname "$SOURCEFEATURES_GROUPED_JSON")"
mkdir -p "$(dirname "$SUBTYPE_ALLSTUDIES_JSON")" "$(dirname "$SUBTYPE_WITH_GT_AA_JSON")"
mkdir -p "$(dirname "$COMPLETEPROFILES_JSON")" "$(dirname "$GT_RAS_JSON")" "$(dirname "$SUBTYPE_RAS_JSON")"

mkdir -p "$OUTPUT_DIR"

cleanup() {
  if [[ -n "${AA_TMP_WORKBOOK:-}" && -f "${AA_TMP_WORKBOOK:-}" ]]; then
    rm -f "$AA_TMP_WORKBOOK"
  fi
}
trap cleanup EXIT

rm -rf "$STAGE_DIR"
mkdir -p "$STAGE_DIR"
rm -f "$REPO_ROOT/temp/build_ns3_sourcefeatures_csv/NS3_SourceFeatures.csv"
rm -f "$REPO_ROOT/temp/build_ns3_sourcefeatures_grouped_csv/NS3_SourceFeatures_Grouped.csv"

"$PYTHON_BIN" "$REPO_ROOT/excel-refid-fasta-discovery/scripts/find_refid_fastas.py" \
  --excel-file "$EXCEL_FILE" \
  --sheet "$SHEET_NAME" \
  --fasta-dir "$FASTA_POOL" \
  --output-dir "$DISCOVERY_TMP" \
  --numpatients-column 'Num Pts' \
  --positive-column NS3Count \
  > "$DISCOVERY_JSON"

DISCOVERY_DIR="$(find "$DISCOVERY_TMP" -maxdepth 1 -type d -name 'refid_fasta_*' | head -n 1)"
cp "$DISCOVERY_DIR/matched_fasta_files.txt" "$MATCHED_TXT"

while IFS= read -r src; do
  [[ -n "$src" ]] || continue
  cp "$src" "$STAGE_DIR/"
done < "$MATCHED_TXT"

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns3_gt_allstudies.py" \
  --excel-file "$EXCEL_FILE" \
  --sheet "$SHEET_NAME" \
  --fasta-dir "$STAGE_DIR" \
  --reference-fasta "$REFERENCE_FASTA" \
  --output-dir "$OUTPUT_DIR" \
  --refid-column RefID \
  --refname-column RefName \
  --numpatients-column 'Num Pts' \
  --positive-column NS3Count \
  > "$GT_ALLSTUDIES_JSON"

if [[ -n "$GENBANK_DIR" ]]; then
  "$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns3_sourcefeatures_csv.py" \
    --matched-fasta-report "$MATCHED_TXT" \
    --genbank-dir "$GENBANK_DIR" \
    > "$SOURCEFEATURES_JSON"

  "$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns3_sourcefeatures_grouped_csv.py" \
    --gt-workbook "$OUTPUT_DIR/NS3_GT_AllStudies.xlsx" \
    --summary-xlsx "$OUTPUT_DIR/NS3_NumSeqs_Naive_1PP_CoversRAS_ByStudy.xlsx" \
    > "$SOURCEFEATURES_GROUPED_JSON"
else
  echo "GENBANK_DIR not provided; skipping NS3 source-feature extraction and grouped summary steps"
fi

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns3_subtype_allstudies_wseqs.py" \
  --combined-workbook "$OUTPUT_DIR/NS3_GT_AllStudies.xlsx" \
  --fasta-dir "$STAGE_DIR" \
  --subtype-json "$SUBTYPE_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$SUBTYPE_ALLSTUDIES_JSON"

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns3_subtype_with_gt_aa.py" \
  --subtype-workbook "$OUTPUT_DIR/NS3_Subtype_AllStudies_WSeqs.xlsx" \
  --fasta-dir "$STAGE_DIR" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$SUBTYPE_WITH_GT_AA_JSON"

AA_TMP_WORKBOOK="$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["output_workbook"])' "$SUBTYPE_WITH_GT_AA_JSON")"

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns3_completeprofiles_tabspergt.py" \
  --input-workbook "$AA_TMP_WORKBOOK" \
  --output-dir "$OUTPUT_DIR" \
  > "$COMPLETEPROFILES_JSON"

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns3_gt_ras_profiles.py" \
  --gt-profile-workbook "$OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$GT_RAS_JSON"

"$PYTHON_BIN" "$REPO_ROOT/scripts/build_ns3_subtype_ras_profiles.py" \
  --subtype-profile-workbook "$OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$SUBTYPE_RAS_JSON"

echo "NS3 pipeline complete"
echo "matched_fasta_report=$MATCHED_TXT"
echo "output_dir=$OUTPUT_DIR"
echo "temp_root=$TEMP_ROOT"
