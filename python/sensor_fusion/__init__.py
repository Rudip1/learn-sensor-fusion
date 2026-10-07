"""sensor_fusion: GNSS, LiDAR and their fusion for vehicle localization and mapping.

The algorithms live in the C++ library (cpp/); this package exposes them through the compiled module
``_core`` and adds data loading and plotting helpers. Point sets are (N, 3) arrays on the Python side.
"""

import sys as _sys

from ._core import cloud, geo, gnss, nmea, sim
from . import data, io, plotting, tables

__version__ = "0.1.0"

# make "import sensor_fusion.nmea" work as well as "from sensor_fusion import nmea"
for _name in ("cloud", "geo", "gnss", "nmea", "sim"):
    _sys.modules[f"{__name__}.{_name}"] = globals()[_name]

__all__ = ["cloud", "data", "geo", "gnss", "io", "nmea", "plotting", "sim", "tables"]
