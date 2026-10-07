"""Turn lists of decoded records into dictionaries of numpy arrays for plotting (no algorithms here)."""

from __future__ import annotations

import numpy as np


def _opt(value) -> float:
    return np.nan if value is None else float(value)


def epochs_to_arrays(epochs) -> dict[str, np.ndarray]:
    """Columns of a list of ``nmea.Epoch``; missing values become NaN.

    Keys: t (seconds since UTC midnight), lat, lon, alt_msl, geoid_sep, fix_quality, sats_used,
    sats_in_view, hdop, pdop, vdop, fix_type, speed, course, rmc_valid.
    """
    cols = {
        "t": [_opt(e.time_of_day_s) for e in epochs],
        "lat": [_opt(e.latitude_deg) for e in epochs],
        "lon": [_opt(e.longitude_deg) for e in epochs],
        "alt_msl": [_opt(e.altitude_msl_m) for e in epochs],
        "geoid_sep": [_opt(e.geoid_separation_m) for e in epochs],
        "fix_quality": [e.fix_quality for e in epochs],
        "sats_used": [e.satellites_used for e in epochs],
        "sats_in_view": [len(e.satellites) if e.satellites else np.nan for e in epochs],
        "hdop": [_opt(e.hdop) for e in epochs],
        "pdop": [_opt(e.pdop) for e in epochs],
        "vdop": [_opt(e.vdop) for e in epochs],
        "fix_type": [e.fix_type for e in epochs],
        "speed": [_opt(e.speed_mps) for e in epochs],
        "course": [_opt(e.course_deg) for e in epochs],
        "rmc_valid": [e.rmc_valid for e in epochs],
    }
    return {k: np.asarray(v) for k, v in cols.items()}


def sky_per_epoch(epochs) -> list[dict[str, np.ndarray]]:
    """Satellites in view for every epoch, carrying the last complete GSV group forward.

    GSV groups are sent less often than GGA/GSA (every few seconds), so each epoch gets the most recent
    sky. Returns one dict per epoch with arrays prn, azimuth_deg, elevation_deg, snr_dbhz (NaN if not
    tracked) and used (bool: PRN listed in that epoch's GSA). Satellites without azimuth/elevation are
    dropped.
    """
    out, last = [], []
    for e in epochs:
        if e.satellites:
            last = [s for s in e.satellites if s.azimuth_deg is not None and s.elevation_deg is not None]
        used = set(e.used_prns)
        out.append(
            {
                "prn": np.array([s.prn for s in last], dtype=int),
                "azimuth_deg": np.array([s.azimuth_deg for s in last], dtype=float),
                "elevation_deg": np.array([s.elevation_deg for s in last], dtype=float),
                "snr_dbhz": np.array([_opt(s.snr_dbhz) for s in last], dtype=float),
                "used": np.array([s.prn in used for s in last], dtype=bool),
            }
        )
    return out
