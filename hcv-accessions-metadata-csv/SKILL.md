---
name: hcv-accessions-metadata-csv
description: Use this skill when the user has RefID-organized HCV FASTA files and a local GenBank flatfile archive and wants to build an accession metadata CSV containing RefID, accession, NS3/NS5A/NS5B gene presence, source-feature qualifiers, and raw structured-comment text.
---

# HCV Accessions Metadata CSV

Use this skill to create the accession metadata CSV consumed by downstream HCV metadata subtype consensus workflows.

## Entry Point

Run:

```bash
python hcv-accessions-metadata-csv/scripts/build_accessions_metadata_csv.py --pipeline-name ns3
```

The script can load `fasta_pool` and `genbank_dir` from `pipeline.local.toml` using `--pipeline-name`.

Explicit inputs:

```bash
python hcv-accessions-metadata-csv/scripts/build_accessions_metadata_csv.py \
  --fasta-dir /path/to/refid_fastas \
  --genbank-dir /path/to/genbank_seq_files \
  --output-csv temp/build_accessions_metadata_csv/Accessions_metadata.csv
```

## Inputs

- `--fasta-dir`: directory of FASTA files whose names begin with `RefID`
- `--genbank-dir`: local GenBank flatfile archive directory containing `.seq` files
- `HCV.fasta` in the repository root for NS3/NS5A/NS5B BLAST gene detection
- optional `pipeline.local.toml` in the repository root for defaults

## Outputs

Default output:

- `temp/build_accessions_metadata_csv/Accessions_metadata.csv`

The CSV includes:

- `RefID`
- `Accession`
- `NS3`, `NS5A`, `NS5B` gene-presence flags
- `StructuredComment`
- `source_*` columns from GenBank source-feature qualifiers

The script prints a run summary with:

- FASTA file count
- accession count
- extracted GenBank record count
- gene-hit accession count
- removed accession count
- final row count
- missing accession count

## Dependency Warning

Before running this skill, warn the user that it requires:

- RefID-organized FASTA files
- local GenBank `.seq` archive files
- BLAST command-line tools: `makeblastdb`, `blastp`, and `blastx`
- root `HCV.fasta`

If any are missing, stop and tell the user which previous setup step is needed.

## Downstream Use

The output CSV is the default metadata input for `$hcv-metadata-subtype-consensus-workflow`.
