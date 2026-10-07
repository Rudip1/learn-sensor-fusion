#!/usr/bin/env python3
"""Execute the solution notebooks top to bottom.

    python tools/run_notebooks.py                 # execute all, fail on the first error
    python tools/run_notebooks.py --inplace       # also store the outputs in the files
    python tools/run_notebooks.py 01 03           # only notebooks whose name starts with these prefixes
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
SOLUTIONS = ROOT / "2_notebooks" / "solutions"


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    inplace = "--inplace" in sys.argv
    paths = sorted(SOLUTIONS.glob("*.ipynb"))
    if args:
        paths = [p for p in paths if any(p.name.startswith(a) for a in args)]
    for path in paths:
        nb = nbformat.read(path, as_version=4)
        start = time.time()
        client = NotebookClient(nb, timeout=900, kernel_name="python3", resources={"metadata": {"path": str(SOLUTIONS)}})
        client.execute()
        print(f"ok  {path.name}  ({time.time() - start:.1f} s)")
        if inplace:
            for cell in nb.cells:
                cell.metadata.pop("execution", None)
            nbformat.write(nb, path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
