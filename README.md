# GenBank Skills

![GenBank Skills logo](assets/genbank-skills-logo.png)

GenBank Skills is a toolkit for working from GenBank accessions to the things researchers usually need next: the linked paper or submission context, clone and quasispecies patterns within a group of samples, and sequence-level checks through reference alignment and downstream quality review.

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

If a PMID is recovered and you want the paper resources, the agent can also resolve the PMCID and print the PMC PDF and package URLs. If you already have one or more paper files locally, give the agent the file paths or a directory and it will convert each supported file to Markdown with a local parser and scan the sentences for accession mentions.

If you want to align one accession, multiple accessions, or a query FASTA file against all entries in a reference FASTA and report the best-matching gene/range, ask the agent with a prompt like:

```text
Use $genbank-reference-alignment to align accession PV289040 against /path/to/reference.fasta and report the best matched gene and aligned range.
```

For multiple accessions:

```text
Use $genbank-reference-alignment to align accessions PV289040, PV289041, and PV289042 against /path/to/reference.fasta and report the best matched gene and aligned ranges in a CSV file.
```

For a query FASTA file:

```text
Use $genbank-reference-alignment to align the sequences in /path/to/query.fasta against /path/to/reference.fasta and report the best matched genes and aligned ranges.
```

If you want to process a cohort of accessions or a GenBank file into a metadata CSV and person-level summary, ask the agent with a prompt like:

```text
Use $genbank-accession-list-metadata to process accessions PV289040, PV289041, and PV289042 into a cohort CSV, identify the patient/person field, count accessions per person, detect likely quasispecies clones, and report the number of persons.
```

For a local GenBank file:

```text
Use $genbank-accession-list-metadata to process /path/to/cohort.gb into a cohort CSV, identify the patient/person field, count accessions per person, detect likely quasispecies clones, and report the number of persons.
```

After `$genbank-accession-list-metadata` finishes, the agent should ask whether you want to continue with `$genbank-gene-split-alignment` for sequence alignment. It should not start alignment automatically.

If you want to extract nucleotide sequences from GenBank records, align them against a multi-gene reference FASTA, and save one aligned nucleotide FASTA per matched gene, ask the agent with a prompt like:

```text
Use $genbank-gene-split-alignment to process accessions PV289040, PV289041, and PV289042 against ./HCV.fasta, align them to all reference genes, and save per-gene nucleotide FASTA files using accession headers.
```

For a local GenBank file:

```text
Use $genbank-gene-split-alignment to process /path/to/cohort.gb against ./HCV.fasta, align all records to the reference genes, and save per-gene nucleotide FASTA files using accession headers.
```

If you started with `$genbank-accession-list-metadata`, the usual next step is:

```text
Yes, continue with $genbank-gene-split-alignment and align this cohort against ./HCV.fasta.
```

If you do not provide the reference FASTA path, the agent should stop and ask which file path to use.
If you provide a gene, it is treated as an optional filter against the FASTA headers.

The agent will run:

```bash
uv run python genbank-single-accession-extractor/scripts/fetch_genbank_accession.py --accession PV289040
```

For alignment, the agent will run:

```bash
uv run python genbank-reference-alignment/scripts/align_accessions_to_reference.py --reference-fasta /path/to/reference.fasta --accession PV289040
```

With an optional gene filter:

```bash
uv run python genbank-reference-alignment/scripts/align_accessions_to_reference.py --reference-fasta /path/to/reference.fasta --gene NS5B --accession PV289040 --accession PV289041 --accession PV289042
```

With a query FASTA file:

```bash
uv run python genbank-reference-alignment/scripts/align_accessions_to_reference.py --reference-fasta /path/to/reference.fasta --query-fasta /path/to/query.fasta
```

For cohort metadata from accessions:

```bash
uv run python genbank-accession-list-metadata/scripts/extract_cohort_metadata.py --accession PV289040 --accession PV289041 --accession PV289042
```

For cohort metadata from a GenBank file:

```bash
uv run python genbank-accession-list-metadata/scripts/extract_cohort_metadata.py --gb-file /path/to/cohort.gb
```

For gene-split nucleotide alignment from accessions:

```bash
uv run python genbank-gene-split-alignment/scripts/split_align_genbank_records.py --reference-fasta ./HCV.fasta --accession PV289040 --accession PV289041 --accession PV289042
```

For gene-split nucleotide alignment from a GenBank file:

```bash
uv run python genbank-gene-split-alignment/scripts/split_align_genbank_records.py --reference-fasta ./HCV.fasta --gb-file /path/to/cohort.gb
```

By default, alignment outputs are written under `outputs/`, not a separate top-level alignment folder.

For cohort-style workflows, the output folder is job-based for safety rather than a fixed name like `cohort_records`. For example:

```text
outputs/cohort_records_20260507T123456Z_ab12cd34/
```

When `$genbank-gene-split-alignment` is run after `$genbank-accession-list-metadata` on the same source cohort, it reuses that cohort folder and writes:

```text
outputs/cohort_records_20260507T123456Z_ab12cd34/
├── metadata.csv
├── person_summary.csv
├── summary.json
├── summary_report.txt
├── gene_alignment_summary.csv
├── gene_alignment_summary.json
├── trees/
└── alignments/
```

The tree step uses FastTree with `-nt -gtr -boot 1000`. If FastTree is not installed on `PATH`, the agent should stop and tell the user to install it before rerunning.

Tree figures are generated from the Newick output when possible, so you get both:

- `trees/<gene>.treefile`
- `trees/<gene>.png`
- `trees/<gene>.svg`

For multiple accessions, the script writes one batch CSV with a simple filename such as:

```text
outputs/alignment_batch_PV289040_PV289041_PV289042.csv
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

## Included Skills

- `genbank-single-accession-extractor/`: Download one GenBank accession, extract FASTA plus organism/source metadata, and enrich references with PMID lookup.
- `genbank-reference-alignment/`: Align one or more GenBank accessions against all reference FASTA entries and report the best matched gene and aligned ranges.
- `genbank-accession-list-metadata/`: Process one cohort of accessions or a GenBank file into metadata CSV plus person/quasispecies summary.
- `genbank-gene-split-alignment/`: Extract nucleotide sequences from GenBank records, align to nucleotide or amino-acid reference genes, and save one codon-preserving nucleotide FASTA per matched gene.
- `hcv-excel-refid-fasta-discovery/`: Filter HCV workbook rows by RefID/patient counts and find matching FASTA files.
- `hcv-gene-genotype-subtype-ref-alignment/`: Build reusable HCV gene genotype/subtype reference alignments and FASTA files.
- `hcv-accessions-metadata-csv/`: Build accession metadata CSVs from RefID FASTA files and local GenBank archives.
- `hcv-ns3-build-workflow/`: Build NS3 genotype, subtype, source-feature, complete-profile, and RAS reports.
- `hcv-ns5a-build-workflow/`: Build NS5A genotype, subtype, source-feature, complete-profile, and RAS reports.
- `hcv-ns5b-build-workflow/`: Build NS5B genotype, subtype, source-feature, complete-profile, and RAS reports.
- `hcv-metadata-subtype-consensus-workflow/`: Build metadata-driven subtype complete profiles, subtype consensus FASTAs, and consensus-to-genotype alignment reports.

The NS3/NS5A/NS5B build workflow wrappers live inside their skill folders. They still load `.env` and `pipeline.local.toml` from the repository root; keep both configuration files in the base folder. Each build workflow skill carries its own `scripts/load_pipeline_defaults.py` helper and passes it the root TOML path explicitly.
The metadata subtype consensus workflow depends on the NS3/NS5A/NS5B build workflow skills and should warn users to complete the required metadata, FASTA, and reference-prep steps before running.

## Repository Layout

```text
.
├── README.md
├── .gitignore
├── genbank-gene-split-alignment/
│   ├── SKILL.md
│   ├── agents/
│   │   └── openai.yaml
│   └── scripts/
│       └── split_align_genbank_records.py
├── genbank-accession-list-metadata/
│   ├── SKILL.md
│   ├── agents/
│   │   └── openai.yaml
│   └── scripts/
│       └── extract_cohort_metadata.py
├── genbank-reference-alignment/
│   ├── SKILL.md
│   ├── agents/
│   │   └── openai.yaml
│   └── scripts/
│       └── align_accessions_to_reference.py
├── genbank-single-accession-extractor/
│   ├── SKILL.md
│   ├── agents/
│   │   └── openai.yaml
│   └── scripts/
│       └── fetch_genbank_accession.py
├── hcv-excel-refid-fasta-discovery/
├── hcv-gene-genotype-subtype-ref-alignment/
├── hcv-accessions-metadata-csv/
│   ├── SKILL.md
│   ├── agents/
│   │   └── openai.yaml
│   └── scripts/
│       └── build_accessions_metadata_csv.py
├── hcv-ns3-build-workflow/
│   ├── SKILL.md
│   ├── NS3_workflow.svg
│   ├── agents/
│   │   └── openai.yaml
│   └── scripts/
│       ├── run_ns3_pipeline.sh
│       ├── load_pipeline_defaults.py
│       └── build_ns3_*.py
├── hcv-ns5a-build-workflow/
│   ├── SKILL.md
│   ├── NS5A_workflow.svg
│   ├── agents/
│   │   └── openai.yaml
│   └── scripts/
│       ├── run_ns5a_pipeline.sh
│       ├── load_pipeline_defaults.py
│       └── build_ns5a_*.py
├── hcv-ns5b-build-workflow/
│   ├── SKILL.md
│   ├── NS5B_workflow.svg
│   ├── agents/
│   │   └── openai.yaml
│   └── scripts/
│       ├── run_ns5b_pipeline.sh
│       ├── load_pipeline_defaults.py
│       └── build_ns5b_*.py
└── hcv-metadata-subtype-consensus-workflow/
    ├── SKILL.md
    ├── agents/
    │   └── openai.yaml
    └── scripts/
        ├── run_metadata_subtype_consensus_workflow.py
        └── export_subtype_consensus_fasta.py
```

## Notes

- `outputs/` is ignored by git.
- FASTA files in this repository are ignored by git.
- Python cache files are ignored by git.
- The extractor keeps the raw GenBank flatfile together with parsed outputs.
- Missing `isolate` is treated as missing data, not as a script error.
- AI fallback cannot run unless `OPENAI_API_KEY` and a model are provided.
