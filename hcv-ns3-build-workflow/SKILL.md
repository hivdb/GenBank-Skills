---
name: hcv-ns3-build-workflow
description: Use this skill when the user wants to run or inspect the HCV NS3 build scripts that create genotype/subtype study workbooks, source-feature summaries, complete profile workbooks, and genotype/subtype RAS profile reports from filtered study FASTA files.
---

# HCV NS3 Build Workflow

Use this skill for the NS3 high-throughput build workflow after study FASTA files have been selected, usually with `$hcv-excel-refid-fasta-discovery`.

## Workflow Chart

See `NS3_workflow.svg` in this skill folder.

## Script Order

1. `scripts/build_ns3_gt_allstudies.py`
2. `scripts/build_ns3_sourcefeatures_csv.py` if GenBank source files are available
3. `scripts/build_ns3_sourcefeatures_grouped_csv.py` if source features were extracted
4. `scripts/build_ns3_subtype_allstudies_wseqs.py`
5. `scripts/build_ns3_subtype_with_gt_aa.py`
6. `scripts/build_ns3_completeprofiles_tabspergt.py`
7. `scripts/build_ns3_gt_ras_profiles.py`
8. `scripts/build_ns3_subtype_ras_profiles.py`

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

## Inputs

- Excel workbook and worksheet containing `RefID`, `RefName`, patient-count, and `NS3Count` fields
- directory of matched/staged NS3 study FASTA files
- `HCV_GT_RefSeqs.fasta`
- `HCV_Subtype_Refs_By_Genome_NA.json`
- `HCV_GT_Refs_By_Gene_AA.json`
- optional GenBank directory for source-feature extraction

## Outputs

The workflow writes NS3 outputs under `outputs/`, including:

- `NS3_GT_AllStudies.xlsx`
- `NS3_SourceFeatures.csv` and grouped source-feature CSV/XLSX outputs when GenBank files are provided
- `NS3_Subtype_AllStudies_WSeqs.xlsx`
- `NS3_Subtype_With_GT_AA.xlsx`
- `NS3_GT_CompleteProfiles_TabsPerGT.xlsx`
- `NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx`
- `NS3_GT_RAS_Profiles.xlsx`
- `NS3_Subtype_RAS_Profiles.xlsx`

## Operating Rules

- Keep NS3 scripts together in this skill folder.
- Use `scripts/run_ns3_pipeline.sh` for complete runs unless the user asks for one specific build step.
- Keep `.env` and `pipeline.local.toml` in the repository root; do not copy them into this skill folder.
- Preserve the order above because later reports consume earlier workbooks.
- If `GENBANK_DIR` is absent, skip only source-feature extraction and grouped source-feature steps.
