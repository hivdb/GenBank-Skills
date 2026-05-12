---
name: hcv-gene-subtype-reference-prep
description: Use this skill when the user wants a one-time preprocessing workflow that builds HCV NS3, NS5A_NTD, and NS5B reference FASTA files from genotype nucleotide references plus subtype genome nucleotide references, then extracts subtype nucleotide and amino-acid reference sequences for those genes with codon-aware alignment.
---

# HCV Gene Subtype Reference Prep

Use this skill when the task is to prepare reusable HCV genotype and subtype reference FASTA files for `NS3`, `NS5A_NTD`, and `NS5B`.

This is a one-time preprocessing step for later high-throughput genotype/subtype calling workflows.

## Workflow

1. Identify the required FASTA/JSON inputs.
   Require:
   - `--gt-gene-na-fasta`
   - `--subtype-genome-na-json`

2. Run the bundled script:

```bash
python3 hcv-gene-subtype-reference-prep/scripts/build_hcv_gene_subtype_refs.py --gt-gene-na-fasta /path/to/HCV_GT_RefSeqs.fasta --subtype-genome-na-json /path/to/HCV_Subtype_Refs_By_Genome_NA.json
```

3. Review the outputs.
   The script:
   - builds separate genotype nucleotide and translated amino-acid reference FASTA files for `NS3`, `NS5A_NTD`, and `NS5B`
   - maps each subtype to its matching genotype-level reference from `HCV_GT_RefSeqs.fasta`
   - extracts per-subtype nucleotide gene sequences from the subtype genome references
   - writes separate translated amino-acid FASTA files for those extracted subtype gene sequences
   - uses codon-aware alignment for subtype extraction and reporting
   - exports an Excel workbook of per-alignment scores and counts
   - writes a plain-text alignment report so alignment quality can be inspected directly

## Output Contract

The script writes one job directory under `outputs/` containing:

- `hcv_gt_gene_refs_ns3_na.fasta`
- `hcv_gt_gene_refs_ns5a_ntd_na.fasta`
- `hcv_gt_gene_refs_ns5b_na.fasta`
- `hcv_gt_gene_refs_ns3_aa.fasta`
- `hcv_gt_gene_refs_ns5a_ntd_aa.fasta`
- `hcv_gt_gene_refs_ns5b_aa.fasta`
- `hcv_subtype_gene_refs_ns3_na.fasta`
- `hcv_subtype_gene_refs_ns5a_ntd_na.fasta`
- `hcv_subtype_gene_refs_ns5b_na.fasta`
- `hcv_subtype_gene_refs_ns3_aa.fasta`
- `hcv_subtype_gene_refs_ns5a_ntd_aa.fasta`
- `hcv_subtype_gene_refs_ns5b_aa.fasta`
- `alignment_scores.xlsx`: per-gene/per-subtype alignment scores and match counts
- `alignment_views.txt`: human-readable codon-aware alignment views with amino-acid match markers
- `summary.json`: machine-readable summary of extracted genes, frames, and source records

## Operating Rules

- Restrict the genotype nucleotide reference set to the eight HCV genotypes for `NS3`, `NS5A_NTD`, and `NS5B`.
- Treat `HCV{genotype}NS5A` in `HCV_GT_RefSeqs.fasta` as full `NS5A`, with `NS5A_NTD` followed by `NS5A_CTD`; use the N-terminal 213 amino acids for `NS5A_NTD`.
- When extracting subtype gene sequences, only compare a subtype genome with a genotype reference from the same genotype.
- Use codon-aware alignment on translated nucleotide frames to locate each gene within the subtype genome sequence.
- Save sequence headers with gene name, genotype, subtype, accession, and source labels.
- Stop with a clear error if a required gene reference cannot be built or a subtype genome cannot be matched for a requested gene.
