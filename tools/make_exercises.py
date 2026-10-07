#!/usr/bin/env python3
"""Derive the reader's notebooks in 2_notebooks/exercises/ from the solved ones in 2_notebooks/solutions/.

In a solution notebook
  - code between the lines "### BEGIN SOLUTION" and "### END SOLUTION" is replaced by a placeholder
    that raises NotImplementedError, keeping the indentation;
  - markdown between "<!-- BEGIN SOLUTION -->" and "<!-- END SOLUTION -->" is replaced by an answer prompt.
All outputs and execution counts are removed.

    python tools/make_exercises.py           # write exercises/
    python tools/make_exercises.py --check   # fail if exercises/ is out of date (used in CI)
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOLUTIONS = ROOT / "2_notebooks" / "solutions"
EXERCISES = ROOT / "2_notebooks" / "exercises"

CODE_RE = re.compile(r"^([ \t]*)### BEGIN SOLUTION[^\n]*\n.*?^[ \t]*### END SOLUTION[^\n]*(\n|$)", re.M | re.S)
MD_RE = re.compile(r"<!-- BEGIN SOLUTION -->.*?<!-- END SOLUTION -->", re.S)


def strip_source(source: str, cell_type: str) -> str:
    if cell_type == "code":
        return CODE_RE.sub(
            lambda m: f"{m.group(1)}# ✏️ your code here\n{m.group(1)}raise NotImplementedError\n", source
        )
    if cell_type == "markdown":
        return MD_RE.sub("✏️ *Your answer here.*", source)
    return source


def make_exercise(nb: dict) -> dict:
    out = json.loads(json.dumps(nb))
    for cell in out["cells"]:
        src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
        src = strip_source(src, cell["cell_type"])
        cell["source"] = src.splitlines(keepends=True)
        if cell["cell_type"] == "code":
            cell["outputs"] = []
            cell["execution_count"] = None
        cell.get("metadata", {}).pop("execution", None)
    return out


def render(nb: dict) -> str:
    return json.dumps(nb, indent=1, ensure_ascii=False) + "\n"


def main() -> int:
    check = "--check" in sys.argv
    EXERCISES.mkdir(parents=True, exist_ok=True)
    stale = []
    for sol in sorted(SOLUTIONS.glob("*.ipynb")):
        text = render(make_exercise(json.loads(sol.read_text())))
        target = EXERCISES / sol.name
        if check:
            if not target.exists() or target.read_text() != text:
                stale.append(target.name)
        else:
            target.write_text(text)
            print(f"wrote {target.relative_to(ROOT)}")
    if stale:
        print("exercise notebooks out of date (run tools/make_exercises.py):", ", ".join(stale))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
