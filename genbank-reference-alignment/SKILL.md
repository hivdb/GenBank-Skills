---
name: genbank-reference-alignment
description: Use this skill when the user provides one accession, multiple accessions, or a query FASTA file and wants each sequence aligned against all entries in a reference FASTA, with the best-matching reference gene reported either in a single-sequence report file or a batch CSV summary.
---

# GenBank Reference Alignment

Use this skill when the task is to align one or more GenBank accessions against all sequences in a reference FASTA and report the best reference/query match range.

## Workflow

1. Identify the query input.
   Accept one accession, a short list of accessions, or a query FASTA file containing one or more sequences.

2. Identify the reference FASTA path and optional gene filter.
   If the user did not provide the reference FASTA path, ask which file path to use.
   If the user provides a gene, use it as a filter against the reference FASTA headers.

3. Run the bundled script:

```bash
uv run python genbank-reference-alignment/scripts/align_accessions_to_reference.py --reference-fasta /path/to/reference.fasta --accession ACCESSION
```

With multiple accessions:

```bash
uv run python genbank-reference-alignment/scripts/align_accessions_to_reference.py --reference-fasta /path/to/reference.fasta --gene GENE --accession ACC1 --accession ACC2
```

With a query FASTA file:

```bash
uv run python genbank-reference-alignment/scripts/align_accessions_to_reference.py --reference-fasta /path/to/reference.fasta --query-fasta /path/to/query.fasta
```

4. Review the outputs:
   - downloaded GenBank flatfiles when accessions are used
   - query FASTA sequences
   - best-matching reference gene record
   - best alignment ranges on reference and query
   - single-sequence text report or batch CSV

## Output Contract

The script writes accession-specific files under `outputs/` by default:

- `record.gb`: downloaded GenBank flatfile for each accession input
- `sequence.fasta`: query sequence FASTA stored per query
- `alignment.json`: structured alignment result for each query

For one query sequence it also writes:

- `alignment_report.txt`: single-sequence readable report

For multiple query sequences it also writes:

- `alignment_batch_*.csv`: one batch-level query-by-query aligned range summary for that request

## Operating Rules

- If the reference FASTA path is missing, ask for it.
- If the gene is provided, it must be found in the reference FASTA headers.
- Ignore FASTA files in this repository for git tracking.
- Default to cached accession downloads unless the user explicitly asks for a full rerun.
