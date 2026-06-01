---
name: genbank-accession-list-metadata
description: Use this skill when the user provides a set of GenBank accessions or a local GenBank file and wants a cohort-level CSV with accession, submission date, and source-feature isolate metadata, plus person-id grouping, quasispecies clone detection, and person-count summary.
---

# GenBank Accession List Metadata

Use this skill when the task is cohort-level GenBank metadata extraction rather than single-record review.

## Workflow

1. Identify the input mode.
   Accept either:
   - one or more `--accession` values
   - one local `--gb-file` path

2. If accessions are provided, download all GenBank flatfiles first.
   If a local GenBank file is provided, read records from that file directly.

3. Extract per-record metadata into a CSV.
   Always include:
   - accession
   - submission date
   - isolate
   - source feature qualifiers

4. Analyze person-level grouping from the extracted metadata.
   Identify the most likely field containing patient or person ID.
   Count accessions per inferred person ID.
   Flag likely quasispecies clones.
   Report the number of persons.

5. Offer the next step.
   After cohort metadata extraction is complete, ask the user whether to continue with `$genbank-gene-split-alignment`.
   Do not start alignment automatically.
   If the user says yes, continue with the extracted GenBank records or cohort `.gb` input and ask for the reference FASTA path if it was not already provided.

## Run

For accessions:

```bash
uv run python genbank-accession-list-metadata/scripts/extract_cohort_metadata.py --accession ACC1 --accession ACC2
```

For a local GenBank file:

```bash
uv run python genbank-accession-list-metadata/scripts/extract_cohort_metadata.py --gb-file /path/to/file.gb
```

## Output Contract

The script writes one cohort directory under `outputs/` containing:

- `metadata.csv`: one row per GenBank record
- `person_summary.csv`: grouped counts by inferred person ID
- `summary.json`: machine-readable cohort summary
- `summary_report.txt`: concise readable report

The cohort directory name is job-based for safety, for example:

- `outputs/cohort_records_20260507T123456Z_ab12cd34/`

If accessions are used, it also stores:

- `records/<accession>.gb`: downloaded GenBank flatfiles

## Operating Rules

- Preserve the raw GenBank records for accession-driven runs.
- Do not guess a patient/person field when there is no evidence; report that it could not be identified.
- Treat clone-like suffixes as supporting evidence for quasispecies, not absolute truth.
- Default to cached downloaded `.gb` files unless the user explicitly asks for a full rerun.
- After reporting cohort findings, explicitly ask whether the user wants to align the sequences with `$genbank-gene-split-alignment`.
