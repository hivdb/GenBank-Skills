#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
from pathlib import Path
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export accession findings to TXT, DOCX, and XLSX."
    )
    parser.add_argument("--accession", required=True, help="Accession directory name under outputs/")
    parser.add_argument("--output-dir", default="outputs", help="Base output directory")
    return parser.parse_args()


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def read_fasta_header(path: Path) -> str:
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            return line
    return ""


def build_report_lines(base_dir: Path) -> list[str]:
    summary = read_json(base_dir / "summary.json")
    organism = read_json(base_dir / "organism.json")
    source_feature = read_json(base_dir / "source_feature.json")
    references = read_json(base_dir / "references.json")
    resolution = read_json(base_dir / "reference_resolution.json")
    fasta_header = read_fasta_header(base_dir / "sequence.fasta")

    ref = references[0] if references else {}
    openai_ref = ref.get("openai_lookup", {})
    qualifiers = source_feature.get("qualifiers", {})

    lines = [
        "GenBank Accession Findings",
        "",
        f"Accession requested: {summary.get('accession_requested')}",
        f"Accession resolved: {summary.get('accession_resolved')}",
        f"Version: {summary.get('version')}",
        f"Definition: {summary.get('definition')}",
        f"Sequence length: {summary.get('sequence_length')}",
        f"FASTA header: {fasta_header}",
        "",
        "Organism",
        f"Scientific name: {summary.get('organism')}",
        f"Common name: {organism.get('common_name')}",
        f"Simple organism name: {source_feature.get('simple_organism_name')}",
        f"Taxonomy ID: {organism.get('taxonomy_id')}",
        "",
        "Source Feature",
        f"Isolate: {source_feature.get('isolate')}",
        f"Host: {qualifiers.get('host')}",
        f"Geo location: {qualifiers.get('geo_loc_name')}",
        f"Collection date: {qualifiers.get('collection_date')}",
        f"Molecule type: {qualifiers.get('mol_type')}",
        f"Note: {qualifiers.get('note')}",
        "",
        "Reference",
        f"Reference: {ref.get('reference')}",
        f"Reference title: {ref.get('title')}",
        f"Reference authors: {ref.get('authors')}",
        f"Reference journal: {ref.get('journal')}",
        f"Matched PMID: {ref.get('pubmed')}",
        f"Matched DOI: {ref.get('openai_doi') or openai_ref.get('doi')}",
        f"Matched paper title: {openai_ref.get('matched_title')}",
        f"Matched paper authors: {openai_ref.get('matched_authors')}",
        f"Matched source URL: {openai_ref.get('source_url')}",
        f"Matching notes: {openai_ref.get('notes')}",
        "",
        "Resolution Steps",
        f"Step 1 status: {resolution.get('step_1_keep_existing_pubmed', [{}])[0].get('status')}",
        f"Step 2 status: {resolution.get('step_2_pubmed_search', [{}])[0].get('status')}",
        f"Step 3 status: {resolution.get('step_3_crossref_lookup', [{}])[0].get('status')}",
        f"Step 4 status: {resolution.get('step_4_openai_web_search_fallback', {}).get('status')}",
        f"OpenAI model: {resolution.get('step_4_openai_web_search_fallback', {}).get('model')}",
    ]
    return lines


def build_rows(lines: list[str]) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for line in lines:
        if not line:
            rows.append(("", ""))
        elif ": " in line:
            key, value = line.split(": ", 1)
            rows.append((key, value))
        else:
            rows.append((line, ""))
    return rows


def write_text(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_docx(path: Path, lines: list[str]) -> None:
    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
</Types>
"""
    rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>
"""
    body = []
    for line in lines:
        text = escape(line)
        body.append(
            f"<w:p><w:r><w:t xml:space=\"preserve\">{text}</w:t></w:r></w:p>"
        )
    document = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:body>
    {''.join(body)}
    <w:sectPr/>
  </w:body>
</w:document>
"""
    with ZipFile(path, "w", compression=ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", rels)
        zf.writestr("word/document.xml", document)


def xlsx_col_name(index: int) -> str:
    result = ""
    while index:
        index, rem = divmod(index - 1, 26)
        result = chr(65 + rem) + result
    return result


def make_inline_cell(ref: str, value: str) -> str:
    return (
        f'<c r="{ref}" t="inlineStr"><is><t xml:space="preserve">'
        f"{escape(value)}</t></is></c>"
    )


def write_xlsx(path: Path, rows: list[tuple[str, str]]) -> None:
    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
  <Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
</Types>
"""
    rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>
"""
    workbook = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
    <sheet name="Findings" sheetId="1" r:id="rId1"/>
  </sheets>
</workbook>
"""
    workbook_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
  <Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>
</Relationships>
"""
    styles = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="1"><font><sz val="11"/><name val="Calibri"/></font></fonts>
  <fills count="1"><fill><patternFill patternType="none"/></fill></fills>
  <borders count="1"><border/></borders>
  <cellStyleXfs count="1"><xf/></cellStyleXfs>
  <cellXfs count="1"><xf xfId="0"/></cellXfs>
  <cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>
"""
    sheet_rows = [
        "<row r=\"1\">"
        + make_inline_cell("A1", "Field")
        + make_inline_cell("B1", "Value")
        + "</row>"
    ]
    for idx, (field, value) in enumerate(rows, start=2):
        sheet_rows.append(
            f"<row r=\"{idx}\">"
            f"{make_inline_cell(f'A{idx}', field)}"
            f"{make_inline_cell(f'B{idx}', value)}"
            "</row>"
        )
    worksheet = f"""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <sheetData>
    {''.join(sheet_rows)}
  </sheetData>
</worksheet>
"""
    with ZipFile(path, "w", compression=ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", rels)
        zf.writestr("xl/workbook.xml", workbook)
        zf.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        zf.writestr("xl/styles.xml", styles)
        zf.writestr("xl/worksheets/sheet1.xml", worksheet)


def main() -> int:
    args = parse_args()
    base_dir = Path(args.output_dir) / args.accession
    lines = build_report_lines(base_dir)
    rows = build_rows(lines)

    write_text(base_dir / "findings_report.txt", lines)
    write_docx(base_dir / "findings_report.docx", lines)
    write_xlsx(base_dir / "findings_report.xlsx", rows)
    print(
        json.dumps(
            {
                "output_dir": str(base_dir),
                "txt": str(base_dir / "findings_report.txt"),
                "docx": str(base_dir / "findings_report.docx"),
                "xlsx": str(base_dir / "findings_report.xlsx"),
            },
            indent=2,
            ensure_ascii=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
