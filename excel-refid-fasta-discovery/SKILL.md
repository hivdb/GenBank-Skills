---
name: excel-refid-fasta-discovery
description: Use this skill when the user has an Excel workbook tab with a RefID column and a NumPatients column and wants the RefIDs where NumPatients is a number greater than 0, then wants matching FASTA files found in a local directory by RefID-prefixed filename.
---

# Excel RefID FASTA Discovery

Use this skill when the task is to read one worksheet from an Excel workbook, keep only rows where `NumPatients > 0`, collect the corresponding `RefID` values, and find FASTA files in a local directory whose filenames begin with those `RefID` values.

## Workflow

1. Identify the required inputs.
   Require:
   - `--excel-file`
   - `--sheet`
   - `--fasta-dir`

2. Run the bundled script:

```bash
uv run python excel-refid-fasta-discovery/scripts/find_refid_fastas.py --excel-file /path/to/workbook.xlsx --sheet TabName --fasta-dir /path/to/fasta_dir
```

3. Review the outputs.
   The script reports:
   - rows where `NumPatients` is numeric and greater than 0
   - the filtered `RefID` list
   - matching FASTA filenames whose basename starts with the `RefID`
   - unmatched `RefID` values

4. Offer the next step.
   After reporting the filenames, ask whether the user wants to continue with a separate skill that reads the individual sequences from those files.
   Do not start sequence reading automatically.

## Output Contract

The script writes one job directory under `outputs/` containing:

- `matching_refids.txt`: one `RefID` per line
- `matched_fasta_files.txt`: one matched FASTA path per line
- `summary.json`: machine-readable summary with matches and unmatched `RefID` values

## Operating Rules

- Treat `NumPatients` as valid only when it can be parsed as a number.
- Keep only rows where the parsed numeric value is strictly greater than `0`.
- Treat blank or missing `RefID` values as unusable and skip them.
- Match FASTA files by filename prefix, using the basename only.
- Search recursively under `--fasta-dir`.
- Default columns are `RefID` and `NumPatients`, but allow overrides when the workbook uses different headers.
