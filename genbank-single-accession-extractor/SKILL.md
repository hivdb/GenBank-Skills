---
name: genbank-single-accession-extractor
description: Use this skill when the user provides one GenBank accession and wants the corresponding GenBank record downloaded, converted to FASTA, and mined for reference metadata plus source-feature qualifiers such as isolate.
---

# GenBank Single Accession Extractor

Use this skill when the task starts from one accession and the goal is to fetch the authoritative GenBank flatfile and derive structured outputs from it.

## Workflow

1. Ask for or identify the accession.
   Accept accession strings with or without version suffixes.

2. Run the bundled script:

```bash
uv run python genbank-single-accession-extractor/scripts/fetch_genbank_accession.py --accession ACCESSION
```

3. Review the generated artifacts and report:
   - downloaded GenBank flatfile
   - FASTA sequence
   - translated organism common name when available
   - extracted and enriched references
   - extracted source feature qualifiers
   - extracted isolate value if present

4. If GenBank references do not contain a PMID, the script should:
   - search PubMed using the extracted citation fields
   - cross-check likely citation matches with Crossref
   - optionally use the OpenAI Responses API with web search for unresolved references
   - for `Direct Submission`, combine accession, authors, scientific organism, and common name in the OpenAI search context

5. If a PMID is available, the script should:
   - resolve the PMCID when the article is in PMC
   - print the PMC PDF URL and the OA package URL for the user
   - if the user supplies one or more local paper files, convert them to Markdown with a local parser and scan the sentences for accession mentions
   - accept repeated paper file paths or a directory of PDFs, DOCX files, and XLSX files

## Output Contract

The script writes an accession-specific output directory containing:

- `record.gb`: downloaded GenBank flatfile
- `sequence.fasta`: FASTA converted from the record
- `organism.json`: scientific and common-name organism metadata
- `references.json`: parsed and enriched `REFERENCE` blocks
- `source_feature.json`: parsed `source` feature qualifiers
- `reference_resolution.json`: PMID lookup diagnostics from PubMed, Crossref, and optional AI fallback
- `pmc_resources.json`: resolved PMCID and PMC download URLs when available
- `paper_markdown/`: Markdown conversion of each supplied paper file when present
- `paper_accession_hits.json`: sentence-level accession match results from the supplied files
- `summary.json`: concise combined summary

## Operating Rules

- Do not invent values that are absent from the record.
- Treat missing `isolate` as a legitimate absence, not an error.
- Preserve the raw downloaded GenBank record alongside parsed outputs.
- If fetch fails, surface the HTTP or accession error clearly.
- Use `uv run python ...` for script execution.
- AI fallback requires `OPENAI_API_KEY` and an explicit model, either `OPENAI_MODEL` or `--openai-model`.
