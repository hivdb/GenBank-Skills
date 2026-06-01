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
2. `scripts/build_ns5a_gt_allstudies.py`
3. `scripts/build_ns5a_sourcefeatures_csv.py` if GenBank source files are available
4. `scripts/build_ns5a_sourcefeatures_grouped_csv.py` if source features were extracted
5. `scripts/build_ns5a_subtype_allstudies_wseqs.py`
6. `scripts/build_ns5a_subtype_with_gt_aa.py`
7. `scripts/build_ns5a_completeprofiles_tabspergt.py`
8. `scripts/build_ns5a_gt_ras_profiles.py`
9. `scripts/build_ns5a_subtype_ras_profiles.py`

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
- optional GenBank directory for source-feature extraction

## Outputs

The workflow writes NS5A outputs under `outputs/`, including:

- `NS5A_GT_AllStudies.xlsx`
- `NS5A_matched_fasta_files.txt`
- `NS5A_SourceFeatures.csv` and grouped source-feature CSV/XLSX outputs when GenBank files are provided
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
- If `GENBANK_DIR` is absent, skip only source-feature extraction and grouped source-feature steps.
