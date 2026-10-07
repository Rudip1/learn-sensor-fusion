"""Readers for the files in data/ (no algorithms here)."""

from __future__ import annotations

import numpy as np

from . import data


def load_navsatfix(relative: str) -> dict[str, np.ndarray]:
    """Read a NavSatFix log converted to CSV (``data/parking_lot/run*_fix.csv``).

    Returns arrays: stamp_ns, t (seconds since the first message), lat, lon, alt (ellipsoidal height),
    status, service, var_e, var_n, var_u (position variances in m^2, east/north/up as in the message).
    """
    raw = np.genfromtxt(data.path(relative), delimiter=",", names=True, dtype=None, encoding="ascii")
    stamp = raw["stamp_ns"].astype(np.int64)
    return {
        "stamp_ns": stamp,
        "t": (stamp - stamp[0]) * 1e-9,
        "lat": raw["latitude_deg"].astype(float),
        "lon": raw["longitude_deg"].astype(float),
        "alt": raw["altitude_m"].astype(float),
        "status": raw["status"].astype(int),
        "service": raw["service"].astype(int),
        "var_e": raw["var_east_m2"].astype(float),
        "var_n": raw["var_north_m2"].astype(float),
        "var_u": raw["var_up_m2"].astype(float),
    }
