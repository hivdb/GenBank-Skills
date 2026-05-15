#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError as exc:  # pragma: no cover
    raise SystemExit(f"tomllib is required: {exc}") from exc


PATH_KEYS = {
    "excel_file",
    "fasta_pool",
    "genbank_dir",
    "output_dir",
    "reference_fasta",
    "subtype_json",
    "gt_aa_json",
    "python_bin",
    "temp_root",
}

KEY_ORDER = [
    "excel_file",
    "fasta_pool",
    "genbank_dir",
    "sheet_name",
    "output_dir",
    "reference_fasta",
    "subtype_json",
    "gt_aa_json",
    "min_sequences",
    "python_bin",
    "temp_root",
]


def shell_escape_double_quoted(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"').replace("$", "\\$").replace("`", "\\`")


def shell_assign(name: str, value: str) -> str:
    return f': "${{{name}:={shell_escape_double_quoted(value)}}}"'


def coerce_value(key: str, value: object, repo_root: Path) -> str:
    text = str(value)
    if key not in PATH_KEYS or not text:
        return text
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = repo_root / path
    return str(path)


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: load_pipeline_defaults.py <pipeline-name> <repo-root>")

    pipeline_name = sys.argv[1]
    repo_root = Path(sys.argv[2]).resolve()
    config_path = repo_root / "pipeline.local.toml"
    if not config_path.is_file():
        return 0

    with config_path.open("rb") as handle:
        data = tomllib.load(handle)

    merged: dict[str, object] = {}
    for section_name in ("common", pipeline_name):
        section = data.get(section_name, {})
        if isinstance(section, dict):
            merged.update(section)

    lines: list[str] = []
    for key in KEY_ORDER:
        value = merged.get(key)
        if value is None:
            continue
        env_name = key.upper()
        lines.append(shell_assign(env_name, coerce_value(key, value, repo_root)))

    if lines:
        print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
