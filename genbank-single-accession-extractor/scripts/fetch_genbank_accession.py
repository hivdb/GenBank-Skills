#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from typing import Any
from dotenv import load_dotenv
from urllib.error import HTTPError, URLError
from urllib.parse import quote_plus, urlencode
from urllib.request import Request, urlopen


EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
CROSSREF_URL = "https://api.crossref.org/works"
OPENAI_RESPONSES_URL = "https://api.openai.com/v1/responses"


load_dotenv()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Download a GenBank record by accession and extract FASTA, "
            "references, and source/isolate feature details."
        )
    )
    parser.add_argument("--accession", required=True, help="GenBank accession or accession.version")
    parser.add_argument(
        "--output-dir",
        default="outputs",
        help="Directory where accession-specific artifacts will be written",
    )
    parser.add_argument(
        "--email",
        default="",
        help="Optional contact email to send to NCBI E-utilities",
    )
    parser.add_argument(
        "--tool",
        default="genbank-single-accession-extractor",
        help="Tool name to send to NCBI E-utilities",
    )
    parser.add_argument(
        "--openai-model",
        default=os.environ.get("OPENAI_MODEL", ""),
        help="Optional OpenAI model for AI PMID fallback lookup",
    )
    parser.add_argument(
        "--disable-ai-lookup",
        action="store_true",
        help="Disable AI-based accession/reference lookup fallback",
    )
    parser.add_argument(
        "--fully-rerun",
        action="store_true",
        help="Ignore existing cached files for this accession and recompute everything",
    )
    return parser.parse_args()


def fetch_genbank_record(accession: str, email: str, tool: str) -> str:
    query = {
        "db": "nuccore",
        "id": accession,
        "rettype": "gbwithparts",
        "retmode": "text",
        "tool": tool,
    }
    if email:
        query["email"] = email
    url = f"{EFETCH_URL}?{urlencode(query)}"

    try:
        with urlopen(url) as response:
            payload = response.read().decode("utf-8")
    except HTTPError as exc:
        raise RuntimeError(f"NCBI request failed with HTTP {exc.code} for accession {accession}") from exc
    except URLError as exc:
        raise RuntimeError(f"NCBI request failed for accession {accession}: {exc.reason}") from exc

    if not payload.strip():
        raise RuntimeError(f"Empty response returned for accession {accession}")
    if "Error occurred" in payload or "Cannot process ID list" in payload:
        raise RuntimeError(f"NCBI returned an accession error for {accession}: {payload.strip()}")

    return payload


def fetch_json(url: str, headers: dict[str, str] | None = None) -> dict[str, Any]:
    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "genbank-single-accession-extractor/1.0",
            **(headers or {}),
        },
    )
    try:
        with urlopen(request) as response:
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        raise RuntimeError(f"Request failed with HTTP {exc.code}: {url}") from exc
    except URLError as exc:
        raise RuntimeError(f"Request failed for {url}: {exc.reason}") from exc


def post_json(url: str, payload: dict[str, Any], headers: dict[str, str] | None = None) -> dict[str, Any]:
    request = Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": "genbank-single-accession-extractor/1.0",
            **(headers or {}),
        },
        method="POST",
    )
    try:
        with urlopen(request) as response:
            return json.loads(response.read().decode("utf-8"))
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"POST failed with HTTP {exc.code}: {detail}") from exc
    except URLError as exc:
        raise RuntimeError(f"POST failed for {url}: {exc.reason}") from exc


def parse_genbank_record(text: str) -> dict[str, Any]:
    lines = text.splitlines()
    metadata = parse_top_level_metadata(lines)
    sequence = parse_origin_sequence(lines)
    references = parse_reference_blocks(lines)
    source_feature = parse_source_feature(lines)
    isolate = source_feature.get("isolate")

    return {
        "metadata": metadata,
        "sequence": sequence,
        "references": references,
        "source_feature": source_feature,
        "isolate": isolate,
    }


def parse_top_level_metadata(lines: list[str]) -> dict[str, str]:
    fields = {
        "locus": "",
        "definition": "",
        "accession": "",
        "version": "",
        "organism": "",
    }
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("FEATURES"):
            break
        if line.startswith("LOCUS"):
            fields["locus"] = line[12:].strip()
        elif line.startswith("DEFINITION"):
            fields["definition"] = collect_continuation(lines, i, 12)
        elif line.startswith("ACCESSION"):
            fields["accession"] = line[12:].strip()
        elif line.startswith("VERSION"):
            fields["version"] = line[12:].strip()
        elif line.startswith("  ORGANISM"):
            fields["organism"] = line[12:].strip()
        i += 1
    return fields


def collect_continuation(lines: list[str], start_index: int, column: int) -> str:
    parts = [lines[start_index][column:].strip()]
    i = start_index + 1
    while i < len(lines):
        next_line = lines[i]
        if len(next_line) < column or next_line[:column].strip():
            break
        parts.append(next_line[column:].strip())
        i += 1
    return " ".join(part for part in parts if part).strip()


def parse_origin_sequence(lines: list[str]) -> str:
    sequence_lines: list[str] = []
    in_origin = False
    for line in lines:
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if not in_origin:
            continue
        if line.startswith("//"):
            break
        sequence_lines.append(re.sub(r"[^A-Za-z]", "", line))
    return "".join(sequence_lines).upper()


def parse_reference_blocks(lines: list[str]) -> list[dict[str, str]]:
    references: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    current_key: str | None = None

    for line in lines:
        if line.startswith("REFERENCE"):
            if current:
                references.append(current)
            current = {"reference": line[12:].strip()}
            current_key = "reference"
            continue

        if current is None:
            continue

        if line.startswith("FEATURES"):
            references.append(current)
            break

        if re.match(r"^(  AUTHORS|  TITLE|  JOURNAL|   PUBMED|  CONSRTM)", line):
            label = line[:12].strip().lower()
            current_key = normalize_reference_key(label)
            current[current_key] = line[12:].strip()
            continue

        if line.startswith("            ") and current_key:
            continuation = line[12:].strip()
            if continuation:
                current[current_key] = f"{current[current_key]} {continuation}".strip()

    if current and (not references or references[-1] is not current):
        references.append(current)

    return references


def normalize_reference_key(label: str) -> str:
    if label == "pubmed":
        return "pubmed"
    return label


def parse_source_feature(lines: list[str]) -> dict[str, Any]:
    feature_start = next((idx for idx, line in enumerate(lines) if line.startswith("FEATURES")), None)
    if feature_start is None:
        return {}

    source_start = None
    for idx in range(feature_start + 1, len(lines)):
        line = lines[idx]
        if line.startswith("ORIGIN"):
            break
        if line.startswith("     source"):
            source_start = idx
            break
    if source_start is None:
        return {}

    qualifiers: dict[str, Any] = {
        "location": lines[source_start][21:].strip(),
        "qualifiers": {},
    }

    current_key: str | None = None
    for idx in range(source_start + 1, len(lines)):
        line = lines[idx]
        if line.startswith("ORIGIN") or re.match(r"^     \S", line):
            break
        qualifier_text = line[21:].rstrip()
        stripped = qualifier_text.strip()
        if not stripped:
            continue

        if stripped.startswith("/"):
            match = re.match(r'^/([^=]+)(?:=(.*))?$', stripped)
            if not match:
                continue
            key = match.group(1)
            value = match.group(2) or ""
            value = clean_qualifier_value(value)
            current_key = key
            add_qualifier_value(qualifiers["qualifiers"], key, value)
            continue

        if current_key:
            continuation = clean_qualifier_value(stripped)
            existing = qualifiers["qualifiers"][current_key]
            if isinstance(existing, list):
                existing[-1] = f"{existing[-1]} {continuation}".strip()
            else:
                qualifiers["qualifiers"][current_key] = f"{existing} {continuation}".strip()

    isolate = qualifiers["qualifiers"].get("isolate")
    if isinstance(isolate, list):
        isolate = isolate[0] if isolate else None
    qualifiers["isolate"] = isolate
    return qualifiers


def clean_qualifier_value(value: str) -> str:
    value = value.strip()
    if value.startswith('"') and value.endswith('"') and len(value) >= 2:
        return value[1:-1]
    return value.strip('"')


def add_qualifier_value(qualifiers: dict[str, Any], key: str, value: str) -> None:
    if key not in qualifiers:
        qualifiers[key] = value
        return
    existing = qualifiers[key]
    if isinstance(existing, list):
        existing.append(value)
        return
    qualifiers[key] = [existing, value]


def write_fasta(path: Path, accession: str, definition: str, sequence: str) -> None:
    header = accession if not definition else f"{accession} {definition}"
    wrapped = "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70))
    path.write_text(f">{header}\n{wrapped}\n", encoding="utf-8")


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def emit_status(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, ensure_ascii=True), file=sys.stderr)


def load_export_module():
    script_path = Path(__file__).with_name("export_accession_report.py")
    spec = spec_from_file_location("export_accession_report", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load export module from {script_path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def export_reports(output_dir: Path, accession: str) -> None:
    module = load_export_module()
    lines = module.build_report_lines(output_dir)
    rows = module.build_rows(lines)
    module.write_text(output_dir / "findings_report.txt", lines)
    module.write_docx(output_dir / "findings_report.docx", lines)
    module.write_xlsx(output_dir / "findings_report.xlsx", rows)
    emit_status(
        {
            "stage": "reports_exported",
            "txt": str(output_dir / "findings_report.txt"),
            "docx": str(output_dir / "findings_report.docx"),
            "xlsx": str(output_dir / "findings_report.xlsx"),
        }
    )


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_cached_run(output_dir: Path) -> dict[str, Any] | None:
    required = {
        "summary": output_dir / "summary.json",
        "organism": output_dir / "organism.json",
        "source_feature": output_dir / "source_feature.json",
        "references": output_dir / "references.json",
        "reference_resolution": output_dir / "reference_resolution.json",
    }
    if not all(path.exists() for path in required.values()):
        return None

    return {
        "summary": read_json(required["summary"]),
        "organism": read_json(required["organism"]),
        "source_feature": read_json(required["source_feature"]),
        "references": read_json(required["references"]),
        "reference_resolution": read_json(required["reference_resolution"]),
    }


def report_files_exist(output_dir: Path) -> bool:
    return all(
        (output_dir / name).exists()
        for name in ("findings_report.txt", "findings_report.docx", "findings_report.xlsx")
    )


def normalize_text(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()


def first_author_surname(authors: str) -> str:
    first = authors.split(" and ")[0].split(",")[0].strip()
    return re.sub(r"[^A-Za-z0-9-]", "", first)


def build_simple_organism_name(organism_info: dict[str, Any]) -> str:
    acronym = str(organism_info.get("acronym") or "").strip()
    if acronym:
        return acronym

    common_name = str(organism_info.get("common_name") or "").strip()
    if common_name:
        if common_name.lower() == "hepatitis c virus":
            return "HCV"
        return common_name

    scientific_name = str(organism_info.get("scientific_name") or "").strip()
    if scientific_name.lower() in {
        "orthohepacivirus hominis",
        "hepacivirus hominis",
        "hepatitis c virus",
    }:
        return "HCV"
    if not scientific_name:
        return ""

    tokens = scientific_name.split()
    if len(tokens) >= 2:
        return " ".join(tokens[:2])
    return scientific_name


def build_pubmed_queries(
    reference: dict[str, str],
    accession: str,
    organism_info: dict[str, Any],
) -> list[str]:
    title = reference.get("title", "").strip()
    authors = reference.get("authors", "").strip()
    journal = reference.get("journal", "").strip()
    organism_name = build_simple_organism_name(organism_info)
    queries: list[str] = []
    surname = first_author_surname(authors)

    if title.lower() == "direct submission":
        if surname and organism_name:
            queries.append(f'{surname}[Author] AND "{organism_name}"[All Fields]')
        elif authors and organism_name:
            queries.append(f'"{authors}"[Author] AND "{organism_name}"[All Fields]')
        return queries

    if title and title.lower() != "direct submission":
        if organism_name:
            queries.append(f'"{title}"[Title] AND "{organism_name}"[All Fields]')
        else:
            queries.append(f'"{title}"[Title]')
        if surname:
            if organism_name:
                queries.append(f'"{title}"[Title] AND {surname}[Author] AND "{organism_name}"[All Fields]')
            else:
                queries.append(f'"{title}"[Title] AND {surname}[Author]')
        if journal:
            if organism_name:
                queries.append(f'"{title}"[Title] AND "{journal}"[Journal] AND "{organism_name}"[All Fields]')
            else:
                queries.append(f'"{title}"[Title] AND "{journal}"[Journal]')
        if organism_name:
            queries.append(f'"{title}"[Title] AND {accession}[All Fields] AND "{organism_name}"[All Fields]')
        else:
            queries.append(f'"{title}"[Title] AND {accession}[All Fields]')
        if surname and journal:
            if organism_name:
                queries.append(
                    f'"{title}"[Title] AND {surname}[Author] AND "{journal}"[Journal] AND "{organism_name}"[All Fields]'
                )
            else:
                queries.append(f'"{title}"[Title] AND {surname}[Author] AND "{journal}"[Journal]')

    if organism_name:
        queries.append(f'{accession}[si] AND "{organism_name}"[All Fields]')
        queries.append(f'{accession} AND "{organism_name}"[All Fields]')
    else:
        queries.append(f"{accession}[si]")
        queries.append(accession)
    return queries


def search_pubmed(
    reference: dict[str, str],
    accession: str,
    organism_info: dict[str, Any],
    email: str,
    tool: str,
) -> dict[str, Any] | None:
    for query in build_pubmed_queries(reference, accession, organism_info):
        params = {
            "db": "pubmed",
            "retmode": "json",
            "retmax": "3",
            "sort": "relevance",
            "term": query,
            "tool": tool,
        }
        if email:
            params["email"] = email
        search = fetch_json(f"{ESEARCH_URL}?{urlencode(params)}")
        idlist = search.get("esearchresult", {}).get("idlist", [])
        if not idlist:
            continue

        summary = fetch_pubmed_summary(idlist[0], email, tool)
        if not summary:
            continue

        summary_title = summary.get("title", "")
        ref_title = reference.get("title", "")
        if ref_title and ref_title.lower() != "direct submission":
            if normalize_text(ref_title) not in normalize_text(summary_title):
                continue

        return {
            "pmid": idlist[0],
            "query": query,
            "summary": summary,
        }
    return None


def fetch_pubmed_summary(pmid: str, email: str, tool: str) -> dict[str, Any] | None:
    params = {
        "db": "pubmed",
        "id": pmid,
        "retmode": "json",
        "tool": tool,
    }
    if email:
        params["email"] = email
    summary = fetch_json(f"{ESUMMARY_URL}?{urlencode(params)}")
    result = summary.get("result", {})
    if pmid not in result:
        return None
    return result[pmid]


def resolve_organism_common_name(organism: str, email: str, tool: str) -> dict[str, Any]:
    if not organism:
        return {
            "scientific_name": organism,
            "common_name": None,
            "source": "missing_organism",
        }

    params = {
        "db": "taxonomy",
        "retmode": "json",
        "retmax": "1",
        "term": f"{organism}[Scientific Name]",
        "tool": tool,
    }
    if email:
        params["email"] = email

    search = fetch_json(f"{ESEARCH_URL}?{urlencode(params)}")
    idlist = search.get("esearchresult", {}).get("idlist", [])
    if not idlist:
        return {
            "scientific_name": organism,
            "common_name": None,
            "source": "taxonomy_not_found",
        }

    summary_params = {
        "db": "taxonomy",
        "id": idlist[0],
        "retmode": "json",
        "tool": tool,
    }
    if email:
        summary_params["email"] = email
    summary = fetch_json(f"{ESUMMARY_URL}?{urlencode(summary_params)}")
    result = summary.get("result", {})
    taxon = result.get(idlist[0], {})

    common_name = (
        taxon.get("commonname")
        or taxon.get("genbankcommonname")
        or taxon.get("blastname")
        or None
    )
    acronym = taxon.get("acronym") or None
    return {
        "taxonomy_id": idlist[0],
        "scientific_name": taxon.get("scientificname") or organism,
        "common_name": common_name,
        "acronym": acronym,
        "blast_name": taxon.get("blastname"),
        "source": "ncbi_taxonomy",
    }


def build_crossref_queries(
    reference: dict[str, str],
    accession: str,
    organism_info: dict[str, Any],
) -> list[str]:
    title = reference.get("title", "").strip()
    if not title or title.lower() == "direct submission":
        return []

    authors = reference.get("authors", "").strip()
    journal = reference.get("journal", "").strip()
    organism_name = build_simple_organism_name(organism_info)
    organism_suffix = f" {organism_name}" if organism_name else ""

    queries = [f"{title}{organism_suffix}"]
    if authors:
        queries.append(f"{title} {authors}{organism_suffix}")
    if journal:
        queries.append(f"{title} {journal}{organism_suffix}")
    queries.append(f"{title} {accession}{organism_suffix}")
    if authors and journal:
        queries.append(f"{title} {authors} {journal}{organism_suffix}")

    seen: set[str] = set()
    deduped: list[str] = []
    for query in queries:
        normalized = query.strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        deduped.append(normalized)
    return deduped


def search_crossref(
    reference: dict[str, str],
    accession: str,
    organism_info: dict[str, Any],
) -> dict[str, Any] | None:
    title = reference.get("title", "").strip()
    if not title or title.lower() == "direct submission":
        return None

    for query in build_crossref_queries(reference, accession, organism_info):
        url = f"{CROSSREF_URL}?rows=3&query.bibliographic={quote_plus(query)}"
        payload = fetch_json(url)
        items = payload.get("message", {}).get("items", [])
        if not items:
            continue

        best_match = None
        for item in items:
            titles = item.get("title") or []
            candidate_title = titles[0] if titles else ""
            if not candidate_title:
                continue
            if normalize_text(title) not in normalize_text(candidate_title):
                continue
            best_match = item
            break

        if best_match is None:
            continue

        return {
            "query": query,
            "title": (best_match.get("title") or [""])[0],
            "doi": best_match.get("DOI"),
            "score": best_match.get("score"),
            "journal": (best_match.get("container-title") or [""])[0],
            "published_print": best_match.get("published-print"),
            "published_online": best_match.get("published-online"),
            "url": best_match.get("URL"),
        }

    return None


def extract_openai_text(response: dict[str, Any]) -> str:
    text_parts: list[str] = []
    for item in response.get("output", []):
        if item.get("type") != "message":
            continue
        for content in item.get("content", []):
            if content.get("type") == "output_text":
                text_parts.append(content.get("text", ""))
    return "\n".join(part for part in text_parts if part).strip()


def extract_openai_sources(response: dict[str, Any]) -> list[dict[str, Any]]:
    sources: list[dict[str, Any]] = []
    for item in response.get("output", []):
        if item.get("type") != "web_search_call":
            continue
        action = item.get("action", {})
        sources.extend(action.get("sources", []))
    return sources


def parse_json_block(text: str) -> Any:
    match = re.search(r"```json\s*(.*?)```", text, flags=re.DOTALL)
    if match:
        return json.loads(match.group(1))

    start = text.find("[")
    end = text.rfind("]")
    if start != -1 and end != -1 and end > start:
        return json.loads(text[start : end + 1])

    start = text.find("{")
    end = text.rfind("}")
    if start != -1 and end != -1 and end > start:
        return json.loads(text[start : end + 1])

    raise ValueError("No JSON payload found in OpenAI response text")


def run_openai_reference_lookup(
    accession: str,
    references: list[dict[str, Any]],
    organism_info: dict[str, Any],
    model: str,
) -> dict[str, Any]:
    api_key = os.environ.get("OPENAI_API_KEY", "").strip()
    if not api_key:
        return {
            "status": "skipped",
            "reason": "OPENAI_API_KEY is not set",
        }
    if not model:
        return {
            "status": "skipped",
            "reason": "OPENAI_MODEL is not set and --openai-model was not provided",
        }

    unresolved = [
        build_openai_reference_payload(ref, accession, organism_info)
        for ref in references
        if not ref.get("pubmed")
    ]
    if not unresolved:
        return {
            "status": "skipped",
            "reason": "All references already contain a PMID",
        }

    prompt = (
        "You are resolving missing citation identifiers for GenBank references. "
        "Search the web and return only JSON as an array. "
        "For each input reference, include reference, pmid, doi, matched_title, matched_authors, source_url, and notes. "
        "If no PMID can be supported, set pmid to null. If no DOI can be supported, set doi to null. "
        "Use the input citation fields to determine which PMID, which DOI, or which paper and author match is best supported. "
        "Do not treat the GenBank accession record itself as the matched paper. "
        "Find a paper, preprint, journal article, or other citation-like publication record if one exists. "
        "If the best evidence only supports the GenBank submission page and not a paper, return pmid as null and doi as null. "
        f"GenBank accession: {accession}. "
        f"Virus scientific name: {organism_info.get('scientific_name')}. "
        f"Virus common name: {organism_info.get('common_name')}. "
        f"References: {json.dumps(unresolved, ensure_ascii=True)}"
    )
    payload = {
        "model": model,
        "reasoning": {"effort": "low"},
        "tools": [{"type": "web_search"}],
        "tool_choice": "auto",
        "include": ["web_search_call.action.sources"],
        "input": prompt,
    }
    response = post_json(
        OPENAI_RESPONSES_URL,
        payload,
        headers={"Authorization": f"Bearer {api_key}"},
    )

    text = extract_openai_text(response)
    try:
        parsed = parse_json_block(text)
    except Exception as exc:
        return {
            "status": "error",
            "reason": f"Failed to parse OpenAI response as JSON: {exc}",
            "raw_text": text,
            "sources": extract_openai_sources(response),
        }

    return {
        "status": "completed",
        "model": model,
        "results": parsed,
        "sources": extract_openai_sources(response),
    }


def build_openai_reference_payload(
    reference: dict[str, Any],
    accession: str,
    organism_info: dict[str, Any],
) -> dict[str, Any]:
    title = str(reference.get("title") or "").strip()
    payload = {
        "reference": reference.get("reference"),
        "authors": reference.get("authors"),
        "title": reference.get("title"),
        "journal": reference.get("journal"),
        "search_context": {
            "accession": accession,
            "authors": reference.get("authors"),
            "virus_scientific_name": organism_info.get("scientific_name"),
            "virus_common_name": organism_info.get("common_name"),
            "submission_note": reference.get("journal"),
        },
    }
    if title and title.lower() != "direct submission":
        payload["search_context"]["title"] = title
    return payload


def apply_step_1_keep_existing_pubmed(reference: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    enriched_ref = dict(reference)
    if reference.get("pubmed"):
        enriched_ref["pmid_source"] = "genbank"
        return enriched_ref, {"status": "kept_existing_pubmed", "pubmed": reference.get("pubmed")}
    return enriched_ref, {"status": "missing_pubmed"}


def apply_step_2_pubmed_search(
    reference: dict[str, Any],
    accession: str,
    organism_info: dict[str, Any],
    email: str,
    tool: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    try:
        pubmed_match = search_pubmed(reference, accession, organism_info, email, tool)
    except Exception as exc:
        return dict(reference), {"status": "error", "error": str(exc)}

    enriched_ref = dict(reference)
    if pubmed_match:
        enriched_ref["pubmed"] = pubmed_match["pmid"]
        enriched_ref["pmid_source"] = "pubmed_esearch"
        enriched_ref["pubmed_match"] = pubmed_match
        return enriched_ref, {"status": "matched", "result": pubmed_match}
    return enriched_ref, {"status": "no_match"}


def apply_step_3_crossref_lookup(
    reference: dict[str, Any],
    accession: str,
    organism_info: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    try:
        crossref_match = search_crossref(reference, accession, organism_info)
    except Exception as exc:
        return dict(reference), {"status": "error", "error": str(exc)}

    enriched_ref = dict(reference)
    if crossref_match:
        enriched_ref["crossref_match"] = crossref_match
        return enriched_ref, {"status": "matched", "result": crossref_match}
    return enriched_ref, {"status": "no_match"}


def apply_step_4_openai_lookup(
    references: list[dict[str, Any]],
    accession: str,
    organism_info: dict[str, Any],
    openai_model: str,
    disable_ai_lookup: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if disable_ai_lookup:
        return references, {"status": "skipped", "reason": "--disable-ai-lookup was set"}

    unresolved = [ref for ref in references if not ref.get("pubmed")]
    if not unresolved:
        return references, {"status": "skipped", "reason": "All references already contain a PMID"}

    try:
        diagnostics = run_openai_reference_lookup(accession, references, organism_info, openai_model)
    except Exception as exc:
        return references, {"status": "error", "reason": str(exc)}

    if diagnostics.get("status") != "completed":
        return references, diagnostics

    results = diagnostics.get("results", [])
    by_reference = {
        item.get("reference"): item
        for item in results
        if isinstance(item, dict) and item.get("reference")
    }
    enriched: list[dict[str, Any]] = []
    for ref in references:
        enriched_ref = dict(ref)
        ai_match = by_reference.get(ref.get("reference"))
        if ai_match:
            enriched_ref["openai_lookup"] = ai_match
            if not enriched_ref.get("pubmed") and ai_match.get("pmid"):
                enriched_ref["pubmed"] = str(ai_match["pmid"])
                enriched_ref["pmid_source"] = "openai_web_search"
            if ai_match.get("doi") and not enriched_ref.get("crossref_match"):
                enriched_ref["openai_doi"] = ai_match.get("doi")
        enriched.append(enriched_ref)
    return enriched, diagnostics


def enrich_references(
    references: list[dict[str, str]],
    accession: str,
    organism_info: dict[str, Any],
    email: str,
    tool: str,
    openai_model: str,
    disable_ai_lookup: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    enriched: list[dict[str, Any]] = []
    diagnostics: dict[str, Any] = {
        "step_1_keep_existing_pubmed": [],
        "step_2_pubmed_search": [],
        "step_3_crossref_lookup": [],
        "step_4_openai_web_search_fallback": {"status": "not_run"},
    }

    for reference in references:
        enriched_ref, step_1 = apply_step_1_keep_existing_pubmed(reference)
        diagnostics["step_1_keep_existing_pubmed"].append(
            {
                "reference": reference.get("reference"),
                "title": reference.get("title"),
                **step_1,
            }
        )
        if not enriched_ref.get("pubmed"):
            enriched_ref, step_2 = apply_step_2_pubmed_search(
                enriched_ref,
                accession,
                organism_info,
                email,
                tool,
            )
            diagnostics["step_2_pubmed_search"].append(
                {
                    "reference": reference.get("reference"),
                    "title": reference.get("title"),
                    **step_2,
                }
            )
        if not enriched_ref.get("crossref_match"):
            enriched_ref, step_3 = apply_step_3_crossref_lookup(
                enriched_ref,
                accession,
                organism_info,
            )
            diagnostics["step_3_crossref_lookup"].append(
                {
                    "reference": reference.get("reference"),
                    "title": reference.get("title"),
                    **step_3,
                }
            )
        enriched.append(enriched_ref)

    enriched, step_4 = apply_step_4_openai_lookup(
        enriched,
        accession,
        organism_info,
        openai_model,
        disable_ai_lookup,
    )
    diagnostics["step_4_openai_web_search_fallback"] = step_4

    return enriched, diagnostics


def main() -> int:
    args = parse_args()
    accession = args.accession.strip()
    if not accession:
        print("Accession must not be empty", file=sys.stderr)
        return 2

    output_dir = Path(args.output_dir) / accession.replace("/", "_")
    output_dir.mkdir(parents=True, exist_ok=True)

    cached_run = None if args.fully_rerun else load_cached_run(output_dir)
    if cached_run is not None:
        emit_status(
            {
                "stage": "using_cached_outputs",
                "output_dir": str(output_dir),
            }
        )
        if not report_files_exist(output_dir):
            export_reports(output_dir, accession)
        print(json.dumps(cached_run["summary"], indent=2, ensure_ascii=True))
        return 0

    record_path = output_dir / "record.gb"
    if record_path.exists() and not args.fully_rerun:
        emit_status(
            {
                "stage": "using_cached_record",
                "record_path": str(record_path),
            }
        )
        record_text = record_path.read_text(encoding="utf-8")
    else:
        record_text = fetch_genbank_record(accession, args.email, args.tool)
        record_path.write_text(record_text, encoding="utf-8")

    parsed = parse_genbank_record(record_text)
    try:
        organism_info = resolve_organism_common_name(
            parsed["metadata"].get("organism", ""),
            args.email,
            args.tool,
        )
    except Exception as exc:
        organism_info = {
            "scientific_name": parsed["metadata"].get("organism", ""),
            "common_name": None,
            "source": "lookup_error",
            "error": str(exc),
        }

    parsed["metadata"]["common_name"] = organism_info.get("common_name")
    parsed["source_feature"]["common_name"] = organism_info.get("common_name")
    simple_organism_name = build_simple_organism_name(organism_info)
    parsed["metadata"]["simple_organism_name"] = simple_organism_name
    parsed["source_feature"]["simple_organism_name"] = simple_organism_name
    emit_status(
        {
            "stage": "organism_parsed",
            "parsed_organism": parsed["metadata"].get("organism"),
            "common_name": organism_info.get("common_name"),
            "simple_organism_name": simple_organism_name,
        }
    )
    enriched_references, reference_resolution = enrich_references(
        parsed["references"],
        accession,
        organism_info,
        args.email,
        args.tool,
        args.openai_model.strip(),
        args.disable_ai_lookup,
    )

    write_fasta(
        output_dir / "sequence.fasta",
        parsed["metadata"].get("version") or parsed["metadata"].get("accession") or accession,
        parsed["metadata"].get("definition", ""),
        parsed["sequence"],
    )
    write_json(output_dir / "references.json", enriched_references)
    write_json(output_dir / "source_feature.json", parsed["source_feature"])
    write_json(output_dir / "organism.json", organism_info)
    write_json(output_dir / "reference_resolution.json", reference_resolution)

    summary = {
        "accession_requested": accession,
        "accession_resolved": parsed["metadata"].get("accession"),
        "version": parsed["metadata"].get("version"),
        "organism": parsed["metadata"].get("organism"),
        "common_name": organism_info.get("common_name"),
        "definition": parsed["metadata"].get("definition"),
        "sequence_length": len(parsed["sequence"]),
        "reference_count": len(enriched_references),
        "reference_pmids_found": sum(1 for ref in enriched_references if ref.get("pubmed")),
        "isolate": parsed["isolate"],
        "output_dir": str(output_dir),
        "openai_lookup_status": reference_resolution.get("step_4_openai_web_search_fallback", {}).get("status"),
    }
    write_json(output_dir / "summary.json", summary)
    export_reports(output_dir, accession)

    print(json.dumps(summary, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
