---
name: genbank-gene-split-alignment
description: Use this skill when the user provides a set of GenBank accessions or a local GenBank file and wants nucleotide sequences extracted, aligned against a multi-gene reference FASTA that may contain nucleotide or amino-acid genes, then split into per-gene aligned nucleotide FASTA files with accession headers.
---

# GenBank Gene Split Alignment

Use this skill when the task is to take many GenBank records, extract their nucleotide sequences, align each one against all genes in a reference FASTA, and save the projected nucleotide alignments into one FASTA file per matched gene.

This skill is the natural follow-up after `$genbank-cohort-metadata`, but it should only start after the agent asks whether the user wants sequence alignment and the user confirms.

## Workflow

1. Identify the input mode.
   Accept either:
   - one or more `--accession` values
   - one local `--gb-file` path

2. Identify the reference FASTA path.
   If the user did not provide the reference FASTA path, ask which file path to use.
   The reference FASTA entries may be nucleotide genes or amino-acid genes.

3. Run the script:

```bash
uv run python genbank-gene-split-alignment/scripts/split_align_genbank_records.py --reference-fasta /path/to/reference.fasta --accession ACC1 --accession ACC2
```

For a local GenBank file:

```bash
uv run python genbank-gene-split-alignment/scripts/split_align_genbank_records.py --reference-fasta /path/to/reference.fasta --gb-file /path/to/file.gb
```

4. Review the outputs:
   - raw GenBank files for accession-driven runs
   - per-record alignment summary
   - per-gene nucleotide alignment FASTA files
   - per-gene FastTree phylogenies and tree figures

## Output Contract

The script writes one cohort directory under `outputs/` containing:

- `gene_alignment_summary.csv`: one row per input record with matched gene and aligned ranges
- `alignments/<gene>.fasta`: one nucleotide alignment FASTA per matched gene
- `gene_alignment_summary.json`: cohort-level alignment summary
- `trees/<gene>.treefile`: FastTree Newick output with 1000 bootstrap replicates
- `trees/<gene>.png` and `trees/<gene>.svg`: tree figures rendered from the Newick tree

If accessions are used, it also stores:

- `records/<accession>.gb`: downloaded GenBank flatfiles

## Operating Rules

- Extract nucleotide sequence from the GenBank records before alignment.
- Align each record against all reference FASTA entries and pick the best-scoring gene match.
- If the reference entry is amino acid, align in amino-acid space and project back to codon-preserving nucleotide alignment.
- Save only nucleotide alignments.
- Use accession as the FASTA header in the per-gene output files.
- Use reference gene positions when projecting aligned query sequences.
- Default to cached downloaded `.gb` files unless the user explicitly asks for a full rerun.
- When invoked after `$genbank-cohort-metadata`, reuse the existing cohort directory if the same source `.gb` file or accession set is detected.
- After writing the aligned FASTA files, run FastTree with `-nt -gtr -boot 1000` for each gene that has at least two sequences.
- If FastTree is not available on `PATH`, stop and ask the user to install it before continuing.
- Generate a tree figure from the FastTree Newick output when rendering libraries are available.
