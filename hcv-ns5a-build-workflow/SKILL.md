---
name: hcv-ns5a-build-workflow
description: Use this skill when the user wants to run or inspect the HCV NS5A build scripts that discover RefID FASTA files, create genotype/subtype study workbooks, source-feature summaries, complete profile workbooks, and genotype/subtype RAS profile reports.
---

# HCV NS5A Build Workflow

Use this skill for the full NS5A high-throughput build workflow. The first step reads the configured Excel worksheet, discovers matching RefID FASTA files, and stages those files for downstream NS5A build steps.

## Workflow Chart

See `NS5A_workflow.svg` in this skill folder.

## Script Order

1. `scripts/find_refid_fastas.py`
2. copy matched FASTA files to `included_refid_fastas/`
3. `scripts/filter_accessions_metadata_by_fasta.py`
4. `scripts/split_refid_metadata_csv.py`
5. `scripts/filter_refid_fastas_by_metadata.py`
6. `scripts/build_ns5a_gt_allstudies.py`
7. `scripts/build_ns5a_sourcefeatures_csv.py` is currently commented out in the wrapper
8. `scripts/build_ns5a_sourcefeatures_grouped_csv.py` is currently commented out in the wrapper
9. `scripts/build_ns5a_subtype_allstudies_wseqs.py`
10. `scripts/build_ns5a_subtype_with_gt_aa.py`
11. `scripts/build_ns5a_completeprofiles_tabspergt.py`
12. `scripts/build_ns5a_gt_ras_profiles.py`
13. `scripts/build_ns5a_subtype_ras_profiles.py`

Prefer the wrapper when running the full workflow:

```bash
EXCEL_FILE=/path/to/HCV_BlastHits.xlsx FASTA_POOL=/path/to/FASTA hcv-ns5a-build-workflow/scripts/run_ns5a_pipeline.sh
```

The shell wrapper is the skill entry point. Do not add a Python entry point unless the orchestration needs cross-platform behavior or richer argument validation; the current Bash wrapper already handles repository defaults, staging, cleanup, and ordered script execution.

Configuration stays in the repository base folder. The wrapper loads:

1. `.env`
2. `pipeline.local.toml`
3. built-in fallbacks

Explicit environment variables provided by the caller take precedence over `pipeline.local.toml`.
The TOML loader is bundled at `scripts/load_pipeline_defaults.py` and is called with the explicit root config path.
Set `sheet_name` in the `[ns5a]` section of `pipeline.local.toml` to choose the input worksheet for discovery and genotype assignment.
Temporary files and step summaries are written under `temp/hcv-ns5a-build-workflow/`.

## Inputs

- Excel workbook and configured worksheet containing `RefID`, `RefName`, patient-count, and `NS5ACount` fields
- FASTA pool directory containing RefID-prefixed FASTA files
- `HCV_GT_RefSeqs.fasta`
- `HCV_Subtype_Refs_By_Genome_NA.json`
- `HCV_GT_Refs_By_Gene_AA.json`
- `Accessions_metadata.csv` for filtering metadata to accessions present in included FASTA files
- optional GenBank directory for source-feature extraction if the commented source-feature steps are re-enabled

The discovery step keeps rows where `RefID` is present and `Num Pts` is not `Exclude`. It does not filter on `NS5ACount` or `Notes`.
After discovery, the wrapper copies all matched RefID FASTA files into `temp/hcv-ns5a-build-workflow/run_ns5a_pipeline/included_refid_fastas/`. Downstream steps that accept `--fasta-dir` use this copied folder, not the original TOML `fasta_pool`.
The metadata filtering step writes `included_accessions_metadata.csv` and reports any FASTA accessions missing from `Accessions_metadata.csv` in `missing_accessions_from_metadata.txt`; both files live in the parent folder of `included_refid_fastas/`.
The per-RefID metadata split step writes CSVs only for RefIDs that have explicit filters under `refid_metadata/`. Current filters: `17` accession is listed in `17.csv`; `29` source_isolate contains `SCRN`; `50` source_isolate contains `week 0`; `85` accession is listed in `85.csv`; `123` source_isolate does not contain `TF`; `142` source_isolate contains `baseline`; `165` accession is listed in `165.csv`; `192` source_isolate contains `day1`; `288` source_isolate contains `pre`; `346` source_isolate contains `baseline/D0`; `535` accession is listed in `535.csv`; `600` source_isolate does not contain `failure`; `661` source_isolation_source equals `plasma`.
The per-RefID FASTA filtering step reads `refid_metadata/RefID_<RefID>_metadata.csv`, keeps only matching `Accession` records in the corresponding copied FASTA file under `included_refid_fastas/`, and prints per-RefID and total before/after record counts.

## Outputs

The workflow writes NS5A outputs under `outputs/`, including:

- `NS5A_GT_AllStudies.xlsx`
- `NS5A_matched_fasta_files.txt`
- discovery `filtered_rows.xlsx` under `temp/hcv-ns5a-build-workflow/.../find_refid_fastas/...`
- copied included RefID FASTA files under `temp/hcv-ns5a-build-workflow/run_ns5a_pipeline/included_refid_fastas/`
- `included_accessions_metadata.csv`
- `missing_accessions_from_metadata.txt`
- `refid_metadata/RefID_<RefID>_metadata.csv`
- filtered copied RefID FASTA files in `included_refid_fastas/` for RefIDs with metadata filters
- source-feature CSV/XLSX outputs only if the commented source-feature steps are re-enabled
- `NS5A_Subtype_AllStudies_WSeqs.xlsx`
- `NS5A_Subtype_With_GT_AA.xlsx`
- `NS5A_GT_CompleteProfiles_TabsPerGT.xlsx`
- `NS5A_Subtype_CompleteProfiles_TabsPerGT.xlsx`
- `NS5A_GT_RAS_Profiles.xlsx`
- `NS5A_Subtype_RAS_Profiles.xlsx`

## Operating Rules

- Keep NS5A scripts together in this skill folder.
- Use `scripts/run_ns5a_pipeline.sh` for complete runs unless the user asks for one specific build step.
- Keep `.env` and `pipeline.local.toml` in the repository root; do not copy them into this skill folder.
- Keep temporary outputs under `temp/hcv-ns5a-build-workflow/` so they do not mix with other skills.
- Preserve the order above because later reports consume earlier workbooks.
- Source-feature extraction and grouped source-feature steps are currently commented out in the wrapper.
