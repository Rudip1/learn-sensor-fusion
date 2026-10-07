"""Locate the small data files shipped in ``data/``.

From a source checkout (``pip install -e .``) the files are read in place. From a plain install (for
example in Colab) they are downloaded once from the repository into ``~/.cache/sensor_fusion``. Set
``SENSOR_FUSION_DATA`` to point at another copy of the ``data/`` folder.
"""

from __future__ import annotations

import os
import urllib.request
from pathlib import Path

RAW_URL = "https://raw.githubusercontent.com/Rudip1/learn-sensor-fusion/main/data/"


def data_dir() -> Path | None:
    """The local ``data/`` folder, or None when running from an installed wheel."""
    env = os.environ.get("SENSOR_FUSION_DATA")
    if env:
        return Path(env)
    candidate = Path(__file__).resolve().parents[2] / "data"
    return candidate if candidate.is_dir() else None


def cache_dir() -> Path:
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "sensor_fusion"


def path(relative: str) -> Path:
    """Absolute path of ``data/<relative>``; downloads the file if there is no local copy."""
    local = data_dir()
    if local is not None and (local / relative).exists():
        return local / relative
    target = cache_dir() / relative
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        tmp = target.with_suffix(target.suffix + ".part")
        urllib.request.urlretrieve(RAW_URL + relative, tmp)
        tmp.replace(target)
    return target
