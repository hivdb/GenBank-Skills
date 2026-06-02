#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
ENV_PATH="$REPO_ROOT/.env"
CONFIG_PATH="$REPO_ROOT/pipeline.local.toml"

if [[ -f "$ENV_PATH" ]]; then
  set -a
  # shellcheck disable=SC1090
  source "$ENV_PATH"
  set +a
fi

if [[ -f "$CONFIG_PATH" ]]; then
  eval "$(python3 "$SCRIPT_DIR/load_pipeline_defaults.py" ns3 "$CONFIG_PATH" "$REPO_ROOT")"
fi

PYTHON_BIN="${PYTHON_BIN:-$REPO_ROOT/.venv/bin/python}"

usage() {
  cat <<'EOF'
Usage:
  EXCEL_FILE=/path/to/HCV_BlastHits.xlsx SHEET_NAME=NS3_PtGT0_Check FASTA_POOL=/path/to/FASTA [GENBANK_DIR=/path/to/genbank_seq_files] hcv-ns3-build-workflow/scripts/run_ns3_pipeline.sh

Optional environment variables:
  SHEET_NAME
  OUTPUT_DIR
  REFERENCE_FASTA
  SUBTYPE_JSON
  GT_AA_JSON
  PYTHON_BIN
  TEMP_ROOT
  ACCESSIONS_METADATA_CSV

Defaults can also be provided in .env and pipeline.local.toml at the repository root.
EOF
}

EXCEL_FILE="${EXCEL_FILE:-}"
FASTA_POOL="${FASTA_POOL:-}"
GENBANK_DIR="${GENBANK_DIR:-}"
SHEET_NAME="${SHEET_NAME:-}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO_ROOT/outputs}"
REFERENCE_FASTA="${REFERENCE_FASTA:-$REPO_ROOT/HCV_GT_RefSeqs.fasta}"
SUBTYPE_JSON="${SUBTYPE_JSON:-$REPO_ROOT/HCV_Subtype_Refs_By_Genome_NA.json}"
GT_AA_JSON="${GT_AA_JSON:-$REPO_ROOT/HCV_GT_Refs_By_Gene_AA.json}"
ACCESSIONS_METADATA_CSV="${ACCESSIONS_METADATA_CSV:-/Users/kaimingtao/Library/CloudStorage/Dropbox/Shared/Vistas/HCVData/Accessions_metadata.csv}"
SKILL_NAME="hcv-ns3-build-workflow"
TEMP_ROOT="${TEMP_ROOT:-$REPO_ROOT/temp/$SKILL_NAME/$(basename "$0" .sh)}"

if [[ -z "$EXCEL_FILE" || -z "$FASTA_POOL" || -z "$SHEET_NAME" ]]; then
  usage
  exit 1
fi

MATCHED_TXT="$OUTPUT_DIR/NS3_matched_fasta_files.txt"
INCLUDED_FASTA_DIR="$TEMP_ROOT/included_refid_fastas"
DISCOVERY_TMP="$TEMP_ROOT/find_refid_fastas"
DISCOVERY_JSON="$DISCOVERY_TMP/discovery_ns3.json"
SKILL_TEMP_ROOT="$REPO_ROOT/temp/$SKILL_NAME"
GT_ALLSTUDIES_JSON="$SKILL_TEMP_ROOT/build_ns3_gt_allstudies/last_run_summary.json"
SOURCEFEATURES_JSON="$SKILL_TEMP_ROOT/build_ns3_sourcefeatures_csv/last_run_summary.json"
SOURCEFEATURES_GROUPED_JSON="$SKILL_TEMP_ROOT/build_ns3_sourcefeatures_grouped_csv/last_run_summary.json"
SUBTYPE_ALLSTUDIES_JSON="$SKILL_TEMP_ROOT/build_ns3_subtype_allstudies_wseqs/last_run_summary.json"
SUBTYPE_WITH_GT_AA_JSON="$SKILL_TEMP_ROOT/build_ns3_subtype_with_gt_aa/last_run_summary.json"
COMPLETEPROFILES_JSON="$SKILL_TEMP_ROOT/build_ns3_completeprofiles_tabspergt/last_run_summary.json"
GT_RAS_JSON="$SKILL_TEMP_ROOT/build_ns3_gt_ras_profiles/last_run_summary.json"
SUBTYPE_RAS_JSON="$SKILL_TEMP_ROOT/build_ns3_subtype_ras_profiles/last_run_summary.json"
REFID_METADATA_DIR="$TEMP_ROOT/refid_metadata"
mkdir -p "$TEMP_ROOT"
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

rm -rf "$INCLUDED_FASTA_DIR"
mkdir -p "$INCLUDED_FASTA_DIR"
rm -rf "$DISCOVERY_TMP"
mkdir -p "$DISCOVERY_TMP"
rm -f "$SKILL_TEMP_ROOT/build_ns3_sourcefeatures_csv/NS3_SourceFeatures.csv"
rm -f "$SKILL_TEMP_ROOT/build_ns3_sourcefeatures_grouped_csv/NS3_SourceFeatures_Grouped.csv"

"$PYTHON_BIN" "$SCRIPT_DIR/find_refid_fastas.py" \
  --excel-file "$EXCEL_FILE" \
  --sheet "$SHEET_NAME" \
  --fasta-dir "$FASTA_POOL" \
  --output-dir "$DISCOVERY_TMP" \
  --numpatients-column 'Num Pts' \
  > "$DISCOVERY_JSON"

DISCOVERY_DIR="$(find "$DISCOVERY_TMP" -maxdepth 1 -type d -name 'refid_fasta_*' | head -n 1)"
cp "$DISCOVERY_DIR/matched_fasta_files.txt" "$MATCHED_TXT"

while IFS= read -r src; do
  [[ -n "$src" ]] || continue
  cp "$src" "$INCLUDED_FASTA_DIR/"
done < "$MATCHED_TXT"

"$PYTHON_BIN" "$SCRIPT_DIR/filter_accessions_metadata_by_fasta.py" \
  --fasta-dir "$INCLUDED_FASTA_DIR" \
  --metadata-csv "$ACCESSIONS_METADATA_CSV" \
  --output-dir "$TEMP_ROOT"

rm -rf "$REFID_METADATA_DIR"
"$PYTHON_BIN" "$SCRIPT_DIR/split_refid_metadata_csv.py" \
  --input-csv "$TEMP_ROOT/included_accessions_metadata.csv" \
  --output-dir "$REFID_METADATA_DIR" \
  > /dev/null

"$PYTHON_BIN" "$SCRIPT_DIR/filter_refid_fastas_by_metadata.py" \
  --metadata-dir "$REFID_METADATA_DIR" \
  --fasta-dir "$INCLUDED_FASTA_DIR"

"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_gt_allstudies.py" \
  --excel-file "$EXCEL_FILE" \
  --sheet "$SHEET_NAME" \
  --fasta-dir "$INCLUDED_FASTA_DIR" \
  --reference-fasta "$REFERENCE_FASTA" \
  --output-dir "$OUTPUT_DIR" \
  --refid-column RefID \
  --refname-column RefName \
  --numpatients-column 'Num Pts' \
  --positive-column NS3Count \
  > "$GT_ALLSTUDIES_JSON"

echo "Skipping NS3 source-feature extraction and grouped summary steps"
# if [[ -n "$GENBANK_DIR" ]]; then
#   "$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_sourcefeatures_csv.py" \
#     --matched-fasta-report "$MATCHED_TXT" \
#     --genbank-dir "$GENBANK_DIR" \
#     > "$SOURCEFEATURES_JSON"
#
#   "$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_sourcefeatures_grouped_csv.py" \
#     --gt-workbook "$OUTPUT_DIR/NS3_GT_AllStudies.xlsx" \
#     --summary-xlsx "$OUTPUT_DIR/NS3_NumSeqs_Naive_1PP_CoversRAS_ByStudy.xlsx" \
#     > "$SOURCEFEATURES_GROUPED_JSON"
# else
#   echo "GENBANK_DIR not provided; skipping NS3 source-feature extraction and grouped summary steps"
# fi

"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_allstudies_wseqs.py" \
  --combined-workbook "$OUTPUT_DIR/NS3_GT_AllStudies.xlsx" \
  --fasta-dir "$INCLUDED_FASTA_DIR" \
  --subtype-json "$SUBTYPE_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$SUBTYPE_ALLSTUDIES_JSON"

"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_with_gt_aa.py" \
  --subtype-workbook "$OUTPUT_DIR/NS3_Subtype_AllStudies_WSeqs.xlsx" \
  --fasta-dir "$INCLUDED_FASTA_DIR" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$SUBTYPE_WITH_GT_AA_JSON"

AA_TMP_WORKBOOK="$("$PYTHON_BIN" -c 'import json,sys; print(json.load(open(sys.argv[1]))["output_workbook"])' "$SUBTYPE_WITH_GT_AA_JSON")"

"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_completeprofiles_tabspergt.py" \
  --input-workbook "$AA_TMP_WORKBOOK" \
  --output-dir "$OUTPUT_DIR" \
  > "$COMPLETEPROFILES_JSON"

"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_gt_ras_profiles.py" \
  --gt-profile-workbook "$OUTPUT_DIR/NS3_GT_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$GT_RAS_JSON"

"$PYTHON_BIN" "$SCRIPT_DIR/build_ns3_subtype_ras_profiles.py" \
  --subtype-profile-workbook "$OUTPUT_DIR/NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx" \
  --gt-aa-json "$GT_AA_JSON" \
  --output-dir "$OUTPUT_DIR" \
  > "$SUBTYPE_RAS_JSON"

echo "NS3 pipeline complete"
echo "matched_fasta_report=$MATCHED_TXT"
echo "included_fasta_dir=$INCLUDED_FASTA_DIR"
echo "output_dir=$OUTPUT_DIR"
echo "temp_root=$TEMP_ROOT"
