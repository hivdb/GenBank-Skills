# Workflow Findings

## Run Commands

### 1. Discover matched FASTA files for `NS3_May11`

```bash
uv run python hcv-ns3-build-workflow/scripts/find_refid_fastas.py \
  --excel-file /path/to/HCV_BlastHits_2026_04_29.xlsx \
  --sheet NS3_May11 \
  --fasta-dir /path/to/FASTA \
  --numpatients-column NumPts \
  --positive-column NS3Count \
  --exclude-known-quasispecies-refids
```

### 2. Assign `NS3` genotype

If the FASTA directory is on Dropbox or another location where temporary BLAST output should not be written, stage the matched FASTA files into a local temp directory first.

```bash
mkdir -p /path/to/local_stage/hcv_fasta_stage_selected
find /path/to/local_stage/hcv_fasta_stage_selected -mindepth 1 -delete
while IFS= read -r f; do
  ln -sf "$f" "/path/to/local_stage/hcv_fasta_stage_selected/$(basename "$f")"
done < outputs/refid_fasta_<workbook_stem>_NS3_May11/matched_fasta_files.txt
```

Then run:

```bash
./.venv/bin/python hcv-ns3-build-workflow/scripts/build_ns3_gt_allstudies.py \
  --excel-file /path/to/HCV_BlastHits_2026_04_29.xlsx \
  --sheet NS3_May11 \
  --fasta-dir /path/to/local_stage/hcv_fasta_stage_selected \
  --reference-fasta HCV_GT_RefSeqs.fasta \
  --numpatients-column NumPts \
  --positive-column NS3Count
```

### 3. Build the combined `NS3` genotype workbook

The genotype script writes `ns3_gt_distance_master.csv`. Convert that CSV into the combined workbook expected by the subtype script:

```bash
./.venv/bin/python - <<'PY'
import csv
from pathlib import Path
from openpyxl import Workbook

base = Path("outputs/<workbook_stem>_NS3_May11_ns3_gt_distance")
csv_path = base / "ns3_gt_distance_master.csv"
xlsx_path = base / "NS3_Alignments_combined.xlsx"

wb = Workbook(write_only=True)
ws = wb.create_sheet("NS3_Alignments_combined")
with csv_path.open(newline="", encoding="utf-8") as fh:
    reader = csv.reader(fh)
    for row in reader:
        ws.append(row)
wb.save(xlsx_path)
print(xlsx_path)
PY
```

### 4. Assign `NS3` subtype

```bash
./.venv/bin/python hcv-ns3-build-workflow/scripts/build_ns3_subtype_allstudies_wseqs.py \
  --combined-workbook outputs/<workbook_stem>_NS3_May11_ns3_gt_distance/NS3_Alignments_combined.xlsx \
  --fasta-dir /path/to/local_stage/hcv_fasta_stage_selected \
  --subtype-json HCV_Subtype_Refs_By_Genome_NA.json
```

### 5. Build NS3 GT resistance profile PNG

```bash
./.venv/bin/python hcv-ns3-build-workflow/scripts/build_ns3_gt_ras_profiles.py \
  --gt-profile-workbook /path/to/NS3_GT_AA_Profiles.xlsx \
  --gt-aa-json HCV_GT_Refs_By_Gene_AA.json
```

This script writes:

- `NS3_GT_Resistance_Profile_Summary.xlsx`
- `NS3_GT_Resistance_Profile_Summary.png`

### 6. Build NS5A GT resistance profile PNG

```bash
./.venv/bin/python hcv-ns5a-build-workflow/scripts/build_ns5a_gt_ras_profiles.py \
  --gt-profile-workbook /path/to/NS5A_GT_AA_Profiles.xlsx \
  --gt-aa-json HCV_GT_Refs_By_Gene_AA.json
```

Default NS5A positions:

- `24,28,30,31,32,58,92,93`

Override them if needed:

```bash
./.venv/bin/python hcv-ns5a-build-workflow/scripts/build_ns5a_gt_ras_profiles.py \
  --gt-profile-workbook /path/to/NS5A_GT_AA_Profiles.xlsx \
  --gt-aa-json HCV_GT_Refs_By_Gene_AA.json \
  --positions 24,28,30,31,32,58,92,93
```

This script writes:

- `NS5A_GT_Resistance_Profile_Summary.xlsx`
- `NS5A_GT_Resistance_Profile_Summary.png`

### 7. Build NS5B GT resistance profile PNG

```bash
./.venv/bin/python hcv-ns5b-build-workflow/scripts/build_ns5b_gt_ras_profiles.py \
  --gt-profile-workbook /path/to/NS5B_GT_AA_Profiles.xlsx \
  --gt-aa-json HCV_GT_Refs_By_Gene_AA.json
```

Default NS5B positions:

- `159,282,316,320,321,414,446,553,554,556,559,561`

Override them if needed:

```bash
./.venv/bin/python hcv-ns5b-build-workflow/scripts/build_ns5b_gt_ras_profiles.py \
  --gt-profile-workbook /path/to/NS5B_GT_AA_Profiles.xlsx \
  --gt-aa-json HCV_GT_Refs_By_Gene_AA.json \
  --positions 159,282,316,320,321,414,446,553,554,556,559,561
```

This script writes:

- `NS5B_GT_Resistance_Profile_Summary.xlsx`
- `NS5B_GT_Resistance_Profile_Summary.png`

### 8. Build HCV gene and subtype reference sets

```bash
./.venv/bin/python hcv-gene-genotype-subtype-ref-alignment/scripts/build_hcv_gene_subtype_refs.py \
  --gt-gene-na-fasta HCV_GT_RefSeqs.fasta \
  --subtype-genome-na-json HCV_Subtype_Refs_By_Genome_NA.json
```

Important detail:

- `frameshift_refinement` is disabled by default in the current reference-prep script so the run finishes reliably.
- use `--enable-frameshift-refinement` only when you explicitly want the slower refinement pass.
- apply study-level `RefID` exclusions during discovery so downstream GT/subtype scripts use the selected input set as-is.
- replace `/path/to/...` placeholders with your actual local paths.

## Step-by-Step Workflows

### NS3

1. Run `hcv-ns3-build-workflow/scripts/find_refid_fastas.py` with:
   - `--sheet NS3_May11`
   - `--numpatients-column NumPts`
   - `--positive-column NS3Count`
2. Stage the matched FASTA files into a local temp directory.
3. Run `hcv-ns3-build-workflow/scripts/build_ns3_gt_allstudies.py`.
4. Convert `ns3_gt_distance_master.csv` to `NS3_Alignments_combined.xlsx`.
5. Run `hcv-ns3-build-workflow/scripts/build_ns3_subtype_allstudies_wseqs.py`.
6. Run `hcv-ns3-build-workflow/scripts/build_ns3_subtype_with_gt_aa.py`.
7. Run `hcv-ns3-build-workflow/scripts/build_ns3_completeprofiles_tabspergt.py`.
8. Run `hcv-ns3-build-workflow/scripts/build_ns3_gt_ras_profiles.py`.
9. Run `hcv-ns3-build-workflow/scripts/build_ns3_subtype_ras_profiles.py`.

Main outputs:

- `NS3_Alignments_combined.xlsx`
- `NS3_Subtype_Alignments_combined.xlsx`
- `NS3_GT_AA_Profiles.xlsx`
- `NS3_Subtype_AA_Profiles.xlsx`
- `NS3_GT_Resistance_Profile_Summary.xlsx`
- `NS3_GT_Resistance_Profile_Summary.png`
- `NS3_Subtype_Resistance_Profile_Summary.xlsx`

### NS5A

Treat `NS5A_NTD` on the AA/profile side as `NS5A` for workflow purposes.

1. Run `hcv-ns5a-build-workflow/scripts/find_refid_fastas.py` with:
   - `--sheet NS5A_PtGT0_Check`
   - `--numpatients-column NumPts`
   - `--positive-column NS5ACount`
2. Stage the matched FASTA files into a local temp directory.
3. Run `hcv-ns5a-build-workflow/scripts/build_ns5a_gt_allstudies.py`.
4. Convert `ns5a_gt_distance_master.csv` to `NS5A_Alignments_combined.xlsx`.
5. Run `hcv-ns5a-build-workflow/scripts/build_ns5a_subtype_allstudies_wseqs.py`.
6. Run `hcv-ns5a-build-workflow/scripts/build_ns5a_subtype_with_gt_aa.py`.
7. Run `hcv-ns5a-build-workflow/scripts/build_ns5a_completeprofiles_tabspergt.py`.
8. Run `hcv-ns5a-build-workflow/scripts/build_ns5a_gt_ras_profiles.py`.
9. Run `hcv-ns5a-build-workflow/scripts/build_ns5a_subtype_ras_profiles.py`.

Main outputs:

- `NS5A_Alignments_combined.xlsx`
- `NS5A_Subtype_Alignments_combined.xlsx`
- `NS5A_GT_AA_Profiles.xlsx`
- `NS5A_Subtype_AA_Profiles.xlsx`
- `NS5A_GT_Resistance_Profile_Summary.xlsx`
- `NS5A_GT_Resistance_Profile_Summary.png`
- `NS5A_Subtype_Resistance_Profile_Summary.xlsx`

### NS5B

1. Run `hcv-ns5b-build-workflow/scripts/find_refid_fastas.py` with:
   - `--sheet NS5B_PtGT0_Check`
   - `--numpatients-column NumPts`
   - `--positive-column NS5BCount`
2. Stage the matched FASTA files into a local temp directory.
3. Run `hcv-ns5b-build-workflow/scripts/build_ns5b_gt_allstudies.py`.
4. Convert `ns5b_gt_distance_master.csv` to `NS5B_Alignments_combined.xlsx`.
5. Run `hcv-ns5b-build-workflow/scripts/build_ns5b_subtype_allstudies_wseqs.py`.
6. Run `hcv-ns5b-build-workflow/scripts/build_ns5b_subtype_with_gt_aa.py`.
7. Run `hcv-ns5b-build-workflow/scripts/build_ns5b_completeprofiles_tabspergt.py`.
8. Run `hcv-ns5b-build-workflow/scripts/build_ns5b_gt_ras_profiles.py`.
9. Run `hcv-ns5b-build-workflow/scripts/build_ns5b_subtype_ras_profiles.py`.

Main outputs:

- `NS5B_Alignments_combined.xlsx`
- `NS5B_Subtype_Alignments_combined.xlsx`
- `NS5B_GT_AA_Profiles.xlsx`
- `NS5B_Subtype_AA_Profiles.xlsx`
- `NS5B_GT_Resistance_Profile_Summary.xlsx`
- `NS5B_GT_Resistance_Profile_Summary.png`
- `NS5B_Subtype_Resistance_Profile_Summary.xlsx`

## What RefID FASTA Discovery Does

The RefID FASTA discovery step does not assign genotype or subtype. It is now bundled as the first script in each NS build workflow skill.

Its job is to:

1. Read one Excel worksheet.
2. Keep rows where the patient count column is numeric and greater than 0.
3. Optionally require one or more additional columns to also be numeric and greater than 0.
4. Collect the corresponding `RefID` values.
5. Find matching FASTA files in a local directory by filename prefix, where the basename starts with the `RefID`.

The scripts are:

- `hcv-ns3-build-workflow/scripts/find_refid_fastas.py`
- `hcv-ns5a-build-workflow/scripts/find_refid_fastas.py`
- `hcv-ns5b-build-workflow/scripts/find_refid_fastas.py`

Important detail:

- The script is not hardcoded to `NS3`.
- It defaults to `RefID` and `NumPatients`.
- Extra filters such as `NS3` must be provided explicitly through `--positive-column`.

## What Actually Assigns Genotype and Subtype

The implemented calling workflow in this repository is an `NS3`-specific downstream process.

### Step 1: Discover relevant FASTA files

Use the local script for the gene workflow being run, for example:

- `hcv-ns3-build-workflow/scripts/find_refid_fastas.py`

Purpose:

- identify which study FASTA files correspond to filtered spreadsheet rows

This step only selects files. It does not classify sequences.

### Step 2: Assign genotype to each accession

Use:

- `hcv-ns3-build-workflow/scripts/build_ns3_gt_allstudies.py`

Purpose:

- read the spreadsheet
- filter studies using patient-count and positive-column rules
- find the FASTA file for each `RefID`
- extract accession sequences from those FASTA files
- align each sequence against `GT1` through `GT8` `NS3` nucleotide references from `HCV_GT_RefSeqs.fasta`
- compute uncorrected nucleotide distance
- assign the best genotype as the closest `GT`

Inputs used by this script:

- Excel workbook and worksheet
- study FASTA directory
- `HCV_GT_RefSeqs.fasta`

Expected behavior from the notes:

- sequences that do not contain `NS3` are skipped
- short alignments can be skipped using a minimum aligned nucleotide threshold
- study-level `RefID` exclusions should already have been applied upstream during discovery/filtering

Outputs:

- one Excel report per study FASTA file
- a combined genotype workbook used by the next step

Related notes:

- `notes/ns3_gt_distance_workflow_2026-05-11.md`

### Step 3: Assign subtype to each accession

Use:

- `hcv-ns3-build-workflow/scripts/build_ns3_subtype_allstudies_wseqs.py`

Purpose:

- read the combined `NS3` genotype workbook
- for each accession, read the assigned `BestGT`
- only compare that accession against subtype references from the same genotype
- use nucleotide alignment to subtype genome references from `HCV_Subtype_Refs_By_Genome_NA.json`
- assign `ClosestSubtype` as the subtype with the lowest uncorrected nucleotide distance
- also report the next closest subtype distance within the same genotype

Inputs used by this script:

- combined genotype workbook such as `NS3_Alignments_combined.xlsx`
- study FASTA directory
- `HCV_Subtype_Refs_By_Genome_NA.json`

Expected behavior from the notes:

- study-level `RefID` exclusions should already have been applied upstream during discovery/filtering
- alignments shorter than 200 nt are skipped by default
- if multiple reference genomes exist for one subtype, the best hit for that subtype is kept

Outputs:

- a combined subtype assignment workbook with fields including:
  - `RefID`
  - `RefName`
  - `AccessionID`
  - `ClosestGT`
  - `ClosestSubtype`
  - `ClosestSubtypeDistance`
  - `NextClosestSubtypeDistance`
  - `AlignedNT`
  - `NextClosestSubtypeAlignedNT`
- this step does not write a separate alignment-view text file

Related notes:

- `notes/ns3_subtype_distance_workflow_2026-05-11.md`

## Overall Pipeline

The concrete workflow I found is:

1. The selected NS build workflow runs its bundled `scripts/find_refid_fastas.py` first to find study FASTA files from filtered spreadsheet rows.
2. `hcv-ns3-build-workflow/scripts/build_ns3_gt_allstudies.py` assigns `NS3` genotype to each accession.
3. `hcv-ns3-build-workflow/scripts/build_ns3_subtype_allstudies_wseqs.py` assigns `NS3` subtype to each accession using genotype-matched subtype references.

## New NS5A Workflow

I added `NS5A` counterparts to the existing `NS3` scripts:

1. `hcv-ns5a-build-workflow/scripts/build_ns5a_gt_allstudies.py`
   - aligns study sequences against `GT1` to `GT8` `NS5A` references from `HCV_GT_RefSeqs.fasta`
   - assigns `BestGT`
   - writes per-study Excel progress files under `NS5A_Alignments.xlsx`
   - writes a master CSV

2. `hcv-ns5a-build-workflow/scripts/build_ns5a_subtype_allstudies_wseqs.py`
   - reads the combined `NS5A` genotype workbook
   - uses each accession's assigned genotype to restrict subtype comparisons
   - aligns against genotype-matched subtype genome references from `HCV_Subtype_Refs_By_Genome_NA.json`
   - assigns `ClosestSubtype` and `NextClosestSubtype`
   - writes `NS5A_Subtype_Alignments_combined.xlsx`

Associated workflow notes:

- `notes/ns5a_gt_distance_workflow_2026-05-13.md`
- `notes/ns5a_subtype_distance_workflow_2026-05-13.md`

## New NS5B Workflow

I added `NS5B` counterparts to the existing `NS3` and `NS5A` scripts:

1. `hcv-ns5b-build-workflow/scripts/build_ns5b_gt_allstudies.py`
   - aligns study sequences against `GT1` to `GT8` `NS5B` references from `HCV_GT_RefSeqs.fasta`
   - assigns `BestGT`
   - writes per-study Excel progress files under `NS5B_Alignments.xlsx`
   - writes a master CSV

2. `hcv-ns5b-build-workflow/scripts/build_ns5b_subtype_allstudies_wseqs.py`
   - reads the combined `NS5B` genotype workbook
   - uses each accession's assigned genotype to restrict subtype comparisons
   - aligns against genotype-matched subtype genome references from `HCV_Subtype_Refs_By_Genome_NA.json`
   - assigns `ClosestSubtype` and `NextClosestSubtype`
   - writes `NS5B_Subtype_Alignments_combined.xlsx`

Associated workflow notes:

- `notes/ns5b_gt_distance_workflow_2026-05-13.md`
- `notes/ns5b_subtype_distance_workflow_2026-05-13.md`

## Scope

This repository does not implement genotype/subtype assignment in the RefID FASTA discovery step.

The explicit classification workflow I found is specific to `NS3`, not a generic all-gene accession typing pipeline.
