---
name: hcv-ns3-build-workflow
description: Use this skill when the user wants to run or inspect the HCV NS3 build scripts that discover RefID FASTA files, create genotype/subtype study workbooks, source-feature summaries, complete profile workbooks, and genotype/subtype RAS profile reports.
---

# HCV NS3 Build Workflow

Use this skill for the full NS3 high-throughput build workflow. The first step reads the configured Excel worksheet, discovers matching RefID FASTA files, and stages those files for downstream NS3 build steps.

## Workflow Chart

See `NS3_workflow.svg` in this skill folder.

## Script Order

1. `scripts/find_refid_fastas.py`
2. copy matched FASTA files to `included_refid_fastas/`
3. `scripts/filter_accessions_metadata_by_fasta.py`
4. `scripts/split_refid_metadata_csv.py`
5. `scripts/filter_refid_fastas_by_metadata.py`
6. `scripts/build_ns3_gt_allstudies.py`
7. `scripts/build_ns3_sourcefeatures_csv.py` is currently commented out in the wrapper
8. `scripts/build_ns3_sourcefeatures_grouped_csv.py` is currently commented out in the wrapper
9. `scripts/build_ns3_subtype_allstudies_wseqs.py`
10. `scripts/build_ns3_subtype_with_gt_aa.py`
11. `scripts/build_ns3_completeprofiles_tabspergt.py`
12. `scripts/export_ns3_consensus_fasta.py`
13. `scripts/build_ns3_gt_ras_profiles.py`
14. `scripts/build_ns3_subtype_ras_profiles.py`

Prefer the wrapper when running the full workflow:

```bash
EXCEL_FILE=/path/to/HCV_BlastHits.xlsx FASTA_POOL=/path/to/FASTA hcv-ns3-build-workflow/scripts/run_ns3_pipeline.sh
```

The shell wrapper is the skill entry point. Do not add a Python entry point unless the orchestration needs cross-platform behavior or richer argument validation; the current Bash wrapper already handles repository defaults, staging, cleanup, and ordered script execution.

Configuration stays in the repository base folder. The wrapper loads:

1. `.env`
2. `pipeline.local.toml`
3. built-in fallbacks

Explicit environment variables provided by the caller take precedence over `pipeline.local.toml`.
The TOML loader is bundled at `scripts/load_pipeline_defaults.py` and is called with the explicit root config path.
Set `sheet_name` in the `[ns3]` section of `pipeline.local.toml` to choose the input worksheet for discovery and genotype assignment.
Temporary files and step summaries are written under `temp/hcv-ns3-build-workflow/`.

## Inputs

- Excel workbook and configured worksheet containing `RefID`, `RefName`, patient-count, and `NS3Count` fields
- FASTA pool directory containing RefID-prefixed FASTA files
- `HCV_GT_RefSeqs.fasta`
- `HCV_Subtype_Refs_By_Genome_NA.json`
- `HCV_GT_Refs_By_Gene_AA.json`
- `Accessions_metadata.csv` for filtering metadata to accessions present in included FASTA files
- optional GenBank directory for source-feature extraction if the commented source-feature steps are re-enabled

The discovery step keeps rows where `RefID` is present and `Num Pts` is not `Exclude`. It does not filter on `NS3Count` or `Notes`.
After discovery, the wrapper copies all matched RefID FASTA files into `temp/hcv-ns3-build-workflow/run_ns3_pipeline/included_refid_fastas/`. Downstream steps that accept `--fasta-dir` must use this copied folder, not the original TOML `fasta_pool`.
The metadata filtering step writes `included_accessions_metadata.csv` and reports any FASTA accessions missing from `Accessions_metadata.csv` in `missing_accessions_from_metadata.txt`; both files live in the parent folder of `included_refid_fastas/`.
The per-RefID metadata split step writes CSVs only for RefIDs that have explicit filters under `refid_metadata/` and prints filter, kept row count, and total row count. Current filters: `30` source_isolate contains `Day1`; `85` accession is listed in `85.csv`; `142` source_isolate contains `baseline`; `192` source_isolate contains `day 1`; `346` source_isolate contains `baseline/D0`; `499` source_isolate contains `HCC`; `600` source_isolate does not contain `failure`; `661` source_isolation_source equals `plasma`; `884` source_isolate contains `Pre-TH`; `943` source_isolate contains `Day 1`; `1356` source_isolate does not contain `IC`; `2008` source_isolate does not contain `chimpanzee`; `2110` source_isolate contains `T0`; `2116` source_collection_date is before 2011; `2138` source_isolate contains `Week 0`; `2150` source_isolate contains `b`; `2168` source_isolate contains `pre`; `2178` source_isolation_source equals `plasma`; `2227` accession is listed in `2227_Nguyen_(2015)_w_metadata_filtered.csv`.
The per-RefID FASTA filtering step reads `refid_metadata/RefID_<RefID>_metadata.csv`, keeps only matching `Accession` records in the corresponding copied FASTA file under `included_refid_fastas/`, and prints per-RefID and total before/after record counts.

## Outputs

The workflow writes NS3 outputs under `outputs/`, including:

- `NS3_GT_AllStudies.xlsx`
- `NS3_matched_fasta_files.txt`
- discovery `filtered_rows.xlsx` under `temp/hcv-ns3-build-workflow/.../find_refid_fastas/...`
- copied included RefID FASTA files under `temp/hcv-ns3-build-workflow/run_ns3_pipeline/included_refid_fastas/`
- `included_accessions_metadata.csv`
- `missing_accessions_from_metadata.txt`
- `refid_metadata/RefID_<RefID>_metadata.csv`
- filtered copied RefID FASTA files in `included_refid_fastas/` for RefIDs with metadata filters
- source-feature CSV/XLSX outputs only if the commented source-feature steps are re-enabled
- `NS3_Subtype_AllStudies_WSeqs.xlsx`
- `NS3_Subtype_With_GT_AA.xlsx`
- `NS3_GT_CompleteProfiles_TabsPerGT.xlsx`
- `NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx`
- `NS3_GT_Consensus.fasta`
- `NS3_Subtype_Consensus.fasta`
- `NS3_GT_RAS_Profiles.xlsx`
- `NS3_Subtype_RAS_Profiles.xlsx`

## Operating Rules

- Keep NS3 scripts together in this skill folder.
- Use `scripts/run_ns3_pipeline.sh` for complete runs unless the user asks for one specific build step.
- Keep `.env` and `pipeline.local.toml` in the repository root; do not copy them into this skill folder.
- Keep temporary outputs under `temp/hcv-ns3-build-workflow/` so they do not mix with other skills.
- Preserve the order above because later reports consume earlier workbooks.
- Source-feature extraction and grouped source-feature steps are currently commented out in the wrapper.
