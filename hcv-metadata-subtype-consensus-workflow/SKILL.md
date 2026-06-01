---
name: hcv-metadata-subtype-consensus-workflow
description: Use this skill when the user has an accession metadata CSV plus RefID-organized FASTA files and wants to build HCV NS3, NS5A, and NS5B subtype complete-profile workbooks, subtype consensus FASTA files, and consensus-to-genotype amino-acid alignment reports.
---

# HCV Metadata Subtype Consensus Workflow

Use this skill for the metadata-driven subtype consensus workflow. It starts from an accession metadata CSV that marks whether each accession has `NS3`, `NS5A`, and/or `NS5B`, plus a FASTA directory whose filenames start with `RefID`.

## Entry Point

Run:

```bash
python hcv-metadata-subtype-consensus-workflow/scripts/run_metadata_subtype_consensus_workflow.py --fasta-dir /path/to/fasta_dir
```

Common inputs:

```bash
--metadata-csv temp/build_accessions_metadata_csv/Accessions_metadata.csv
--output-root temp/metadata_subtype_consensus_workflow
--subtype-json HCV_Subtype_Refs_By_Genome_NA.json
--gt-aa-json HCV_GT_Refs_By_Gene_AA.json
--gt-aa-fasta HCV_GT_Refs_NS3_NS5A_NTD_NS5B_AA.fasta
--reference-fasta HCV_GT_RefSeqs.fasta
--genes NS3 NS5A NS5B
```

To only rebuild consensus-to-genotype alignment text reports from existing consensus FASTAs:

```bash
python hcv-metadata-subtype-consensus-workflow/scripts/run_metadata_subtype_consensus_workflow.py --only-consensus-alignments
```

## Dependency Warning

Before running this skill, warn the user that it depends on outputs and scripts from earlier HCV workflow steps.

Required previous setup:

- accession metadata CSV exists, usually from `$hcv-accessions-metadata-csv`
- RefID-organized FASTA files are available
- genotype/subtype reference files exist:
  - `HCV_GT_RefSeqs.fasta`
  - `HCV_Subtype_Refs_By_Genome_NA.json`
  - `HCV_GT_Refs_By_Gene_AA.json`
  - `HCV_GT_Refs_NS3_NS5A_NTD_NS5B_AA.fasta`

Required sibling workflow skills:

- `$hcv-ns3-build-workflow`
- `$hcv-ns5a-build-workflow`
- `$hcv-ns5b-build-workflow`

If any of these are missing, stop and tell the user which previous step to run first. The user only needs to know that this skill uses the gene workflow skills internally; do not explain internal script path wiring unless troubleshooting.

## Workflow

For each selected gene:

1. Read `metadata.csv` and select accessions where the gene column is present.
2. Stage RefID FASTA files containing the selected accessions.
3. Create a seed workbook for genotype assignment.
4. Run the gene workflow's genotype assignment script.
5. Run the gene workflow's subtype assignment script.
6. Run the gene workflow's subtype-with-GT-AA extraction script.
7. Run the gene workflow's complete-profile script.
8. Export subtype consensus FASTA with the bundled helper.
9. Write subtype-consensus-to-genotype-AA alignment text reports.

## Outputs

The script writes one `output-root` containing:

- `workflow_summary.json`
- one subdirectory per gene
- `<gene>_Metadata_GT_Seed.xlsx`
- `<gene>_gt_assignment_summary.json`
- `<gene>_subtype_assignment_summary.json`
- `<gene>_aa_extraction_summary.json`
- `<gene>_completeprofiles_summary.json`
- `<gene>_Subtype_Consensus.fasta`
- `<gene>_Subtype_Consensus_GT_Alignment.txt`

## Operating Rules

- Keep `export_subtype_consensus_fasta.py` bundled in this skill because it is only used by this workflow.
- Before running the full workflow, warn the user that `$hcv-ns3-build-workflow`, `$hcv-ns5a-build-workflow`, and `$hcv-ns5b-build-workflow` are required.
- If the user asks to modify genotype, subtype, AA extraction, or complete-profile behavior, edit the owning gene workflow skill, not this orchestrator.
- If the user asks to modify consensus FASTA export or consensus-to-GT alignment reporting, edit this skill.
