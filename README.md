# GenBank Skills

This repository helps an AI agent work with a single GenBank accession at a time.

## Start With Codex

This workflow is intended to be started in Codex.

Recommended starting point:

- Open Codex: `https://chatgpt.com/codex`

Official OpenAI docs:

- Codex overview: `https://platform.openai.com/docs/codex/overview`
- Codex getting started: `https://openai.com/academy/codex-how-to-start`
- Working with Codex: `https://openai.com/academy/working-with-codex/`

## How To Use

If you have a GenBank accession and want the GenBank record, FASTA, references, organism common name, and isolate metadata, ask the agent with a prompt like:

```text
Use $genbank-single-accession-extractor to download accession PV289040, extract FASTA, organism common name, and isolate metadata, and recover PMID information for references when possible.
```

The agent will run:

```bash
uv run python genbank-single-accession-extractor/scripts/fetch_genbank_accession.py --accession PV289040
```

By default, the script reuses existing per-accession intermediate/output files if they already exist under `outputs/ACCESSION/`.

Only use a full recomputation when you explicitly want to ignore cached files:

```bash
uv run python genbank-single-accession-extractor/scripts/fetch_genbank_accession.py --accession PV289040 --fully-rerun
```

This creates:

```text
outputs/PV289040/
├── record.gb
├── sequence.fasta
├── organism.json
├── references.json
├── reference_resolution.json
├── source_feature.json
└── summary.json
```

If a GenBank reference already contains a PMID, it is kept.

If a reference does not contain a PMID, the script will try to resolve it in this order:

1. PubMed search using the extracted title, authors, and accession
2. Crossref citation lookup for likely DOI and title matches
3. Optional OpenAI API web-search fallback for unresolved references

## PMID And DOI Search Strategy

The extractor uses the information already present in the GenBank record and applies a staged matching strategy.

### Step 1: Keep Existing PubMed IDs

If a `REFERENCE` block already includes a `PUBMED` field, that PMID is preserved as the authoritative match.

### Step 2: PubMed Search

For references without a PMID, the script searches PubMed first.

The simple organism name is always included in PubMed Step 2 queries.

For normal article-like references, it builds queries from:

- reference title
- first author surname
- journal
- simple organism name parsed from the GenBank record
- accession

Typical PubMed query patterns are:

- exact title plus simple organism name
- exact title plus first author plus simple organism name
- exact title plus journal plus simple organism name
- exact title plus accession plus simple organism name
- exact title plus first author plus journal plus simple organism name
- accession plus simple organism name

The script then checks the returned PubMed summary and only accepts a match when the title is consistent with the GenBank reference.

If the title is `Direct Submission`, PubMed Step 2 does not use the title or accession query patterns. It only uses:

- first author surname plus simple organism name

If the surname is not usable, it falls back to:

- full author string plus simple organism name

Before reference search begins, the script reports the organism parsed from the GenBank file, including:

- scientific organism name
- common name when available
- simple organism name used in search queries

### Step 3: Crossref Lookup For DOI Context

If PMID is still missing, the script queries Crossref using the same extracted context as PubMed.

The simple organism name is always included in Crossref Step 3 queries.

It builds Crossref lookup strings from:

- title
- first author surname or full author string
- journal
- simple organism name parsed from the GenBank record
- accession

Typical Crossref query patterns are:

- title plus simple organism name
- title plus authors plus simple organism name
- title plus journal plus simple organism name
- title plus accession plus simple organism name
- title plus authors plus journal plus simple organism name

Crossref is used to collect likely DOI and citation candidates. This helps confirm whether the reference behaves like a published article or just a submission note.

### Step 4: OpenAI Web Search Fallback

If PMID is still unresolved and OpenAI API access is configured, the script sends the unresolved references to the OpenAI Responses API with web search enabled.

The OpenAI search context always includes:

- accession
- authors
- organism scientific name
- organism common name when available
- submission note from the GenBank journal field

If the title is not `Direct Submission`, the title is also included in the search context.

The model returns structured JSON with:

- reference
- pmid
- doi
- matched title
- matched authors
- source URL
- notes

The current prompt template used by the script is:

```text
You are resolving missing citation identifiers for GenBank references. Search the web and return only JSON as an array. For each input reference, include reference, pmid, doi, matched_title, matched_authors, source_url, and notes. If no PMID can be supported, set pmid to null. If no DOI can be supported, set doi to null. Use the input citation fields to determine which PMID, which DOI, or which paper and author match is best supported. Do not treat the GenBank accession record itself as the matched paper. Find a paper, preprint, journal article, or other citation-like publication record if one exists. If the best evidence only supports the GenBank submission page and not a paper, return pmid as null and doi as null. GenBank accession: {accession}. Virus scientific name: {scientific_name}. Virus common name: {common_name}. References: {unresolved_references_json}
```

### Matching Outcome

The final decision is conservative:

- if a PMID is supported, it is added to `references.json`
- if a DOI is supported by the OpenAI fallback, it is attached to the reference as `openai_doi`
- if no PMID can be supported, the script leaves it empty instead of guessing
- diagnostics from PubMed, Crossref, and OpenAI are written to `reference_resolution.json`

## What You Get

The agent returns:

- accession and version
- organism and record definition
- translated organism common name when available from NCBI taxonomy
- FASTA output file
- extracted and enriched `REFERENCE` entries
- extracted `source` feature qualifiers
- extracted `isolate` value when present
- PMID recovery status for each reference

## Direct Script Use

You can run the extractor directly:

```bash
uv run python genbank-single-accession-extractor/scripts/fetch_genbank_accession.py --accession ACCESSION
```

Optional flags:

```bash
--output-dir outputs
--email you@example.org
--tool genbank-single-accession-extractor
--openai-model YOUR_MODEL
--disable-ai-lookup
--fully-rerun
```

If you want AI-based accession/reference lookup, set:

```bash
export OPENAI_API_KEY=...
export OPENAI_MODEL=...
```

Without those settings, the script still runs and still does GenBank, PubMed, and Crossref lookup. It only skips the OpenAI fallback.

Cache behavior:

- default: reuse existing cached files for the accession if the expected output JSON files already exist
- default fallback: if `record.gb` already exists but final JSON outputs do not, reuse `record.gb` and recompute the downstream parsing and lookup steps
- `--fully-rerun`: ignore existing files and fetch/recompute everything again

## Included Skill

- `genbank-single-accession-extractor/`: Download one GenBank accession, extract FASTA plus organism/source metadata, and enrich references with PMID lookup.

## Repository Layout

```text
.
├── README.md
├── .gitignore
└── genbank-single-accession-extractor/
│   ├── SKILL.md
│   ├── agents/
│   │   └── openai.yaml
│   └── scripts/
│       └── fetch_genbank_accession.py
```

## Notes

- `outputs/` is ignored by git.
- Python cache files are ignored by git.
- The extractor keeps the raw GenBank flatfile together with parsed outputs.
- Missing `isolate` is treated as missing data, not as a script error.
- AI fallback cannot run unless `OPENAI_API_KEY` and a model are provided.
