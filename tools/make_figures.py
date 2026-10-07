#!/usr/bin/env python3
"""Regenerate the figures in 1_theory/figures/ from the sensor_fusion module.

    python tools/make_figures.py            # all chapters
    python tools/make_figures.py 1 3        # selected chapters

Each chapter function also prints the numbers quoted in the theory text, so they can be checked.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import sensor_fusion as sf
from sensor_fusion import geo, gnss, nmea
from sensor_fusion.plotting import SERIES, figure

OUT = Path(__file__).resolve().parents[1] / "1_theory" / "figures"


def save(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name)
    print(f"  wrote {OUT / name}")


def spherical_receiver(lat_deg: float, lon_deg: float, radius: float = 6_371_000.0) -> np.ndarray:
    lat, lon = np.radians([lat_deg, lon_deg])
    return radius * np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def chapter1() -> None:
    print("chapter 1")
    r = spherical_receiver(47.47, 19.06)
    az = np.radians([20, 70, 115, 160, 210, 250, 300, 340])
    el = np.radians([65, 30, 45, 15, 55, 25, 40, 12])
    sats = gnss.place_satellites(r, az, el)
    bias = 1e-3 * gnss.SPEED_OF_LIGHT
    rho_true = gnss.geometric_ranges(r, sats) + bias
    rng = np.random.default_rng(1)
    rho_noisy = rho_true + rng.normal(0.0, 3.0, len(rho_true))

    fig, ax = figure(width=6.0, height=3.2)
    for rho, label, color in [(rho_true, "noiseless", SERIES[0]), (rho_noisy, r"$\sigma_\rho$ = 3 m", SERIES[1])]:
        errors = [np.linalg.norm(r)]  # start at the Earth's centre
        for k in range(1, 8):
            sol = gnss.solve_position(sats, rho, max_iterations=k, tolerance_m=0.0)
            errors.append(np.linalg.norm(sol.position - r))
        ax.semilogy(range(len(errors)), np.maximum(errors, 1e-9), "o-", color=color, label=label)
        print(f"  {label}: error per iteration [m] = " + ", ".join(f"{e:.3g}" for e in errors))
    ax.set_xlabel("Gauss-Newton iteration")
    ax.set_ylabel("position error [m]")
    ax.set_title("Single-point positioning, 8 satellites, start at the Earth's centre")
    ax.legend()
    save(fig, "01_gauss_newton.png")

    fig, axes = figure(2, 1, width=7.0, height=4.4, sharex=False)
    for ax, name in zip(axes, ["nmea_log1.nmea", "nmea_log2.nmea"]):
        log = nmea.parse_log_file(str(sf.data.path("gnss/" + name)))
        a = sf.tables.epochs_to_arrays(log.epochs)
        t = a["t"] - a["t"][0]
        ax.step(t, a["sats_used"], where="post", color=SERIES[0], label="satellites used (GGA)")
        ax.step(t, a["fix_quality"], where="post", color=SERIES[1], label="fix quality (GGA)")
        ax.step(t, a["fix_type"], where="post", color=SERIES[2], label="fix type (GSA)")
        ax.set_title(name)
        ax.set_ylabel("count / code")
        valid = a["fix_quality"] > 0
        print(
            f"  {name}: epochs {len(t)}, first fix after {t[valid][0]:.0f} s, "
            f"satellites used {a['sats_used'].min()}..{a['sats_used'].max()}"
        )
    axes[1].set_xlabel("time since first epoch [s]")
    axes[0].legend(ncol=3, loc="lower right")
    save(fig, "01_nmea_logs.png")


def chapter2() -> None:
    print("chapter 2")
    # 2.1 meridian section with exaggerated flattening: geodetic vs geocentric latitude, N(phi)
    a, f_vis = 1.0, 0.3
    b = a * (1 - f_vis)
    e2 = f_vis * (2 - f_vis)
    t = np.linspace(0, np.pi / 2, 200)
    fig, ax = figure(width=5.2, height=4.4)
    ax.plot(a * np.cos(t), b * np.sin(t), color=SERIES[0], lw=2)
    lat = np.radians(45.0)
    N = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    px, pz = N * np.cos(lat), N * (1 - e2) * np.sin(lat)
    ax.plot([px, px - N * np.cos(lat)], [pz, pz - N * np.sin(lat)], color=SERIES[1], lw=1.4)
    ax.plot([0, px], [0, pz], color=SERIES[2], lw=1.4)
    ax.plot([px, px + 0.35 * np.cos(lat)], [pz, pz + 0.35 * np.sin(lat)], color=SERIES[1], lw=1.4, ls="--")
    ax.plot(px, pz, "o", color="black", ms=5)
    ax.text(px + 0.12, pz + 0.22, r"ellipsoid normal: geodetic latitude $\varphi$", color=SERIES[1], ha="right")
    ax.text(0.5 * px - 0.05, 0.5 * pz + 0.05, "geocentric radius", color=SERIES[2], ha="right", rotation=26)
    ax.text(0.5 * (px - N * np.cos(lat)) + 0.47 * px, 0.5 * (2 * pz - N * np.sin(lat)) - 0.12,
            r"$N(\varphi)$: surface to polar axis", color=SERIES[1])
    ax.axhline(0, color="0.6", lw=0.8)
    ax.axvline(0, color="0.6", lw=0.8)
    ax.text(1.0, 0.03, "a", color="0.3")
    ax.text(0.03, b + 0.02, "b", color="0.3")
    ax.set_aspect("equal")
    ax.set_xlim(-0.25, 1.25)
    ax.set_ylim(-0.4, 1.0)
    ax.set_xlabel(r"distance from the polar axis $\sqrt{x^2+y^2}$")
    ax.set_ylabel("z")
    ax.set_title("Meridian section (flattening exaggerated to f = 0.3)")
    save(fig, "02_ellipsoid.png")

    # 2.2 error of local approximations vs distance
    lat0, lon0, h0 = 47.4724, 19.0633, 148.0
    ltp = geo.LocalTangentPlane(lat0, lon0, h0)
    d = np.logspace(1, 5, 61)  # 10 m .. 100 km, half north half east
    enu_true = np.column_stack([d / np.sqrt(2), d / np.sqrt(2), np.zeros_like(d)])
    llh = ltp.enu_to_geodetic(enu_true)
    exact = ltp.geodetic_to_enu(llh[:, 0], llh[:, 1], llh[:, 2])
    flat = geo.geodetic_to_enu_flat(llh[:, 0], llh[:, 1], llh[:, 2], lat0, lon0, h0)
    k = 111_000.0
    crude = np.column_stack([(llh[:, 1] - lon0) * k * np.cos(np.radians(lat0)), (llh[:, 0] - lat0) * k,
                             llh[:, 2] - h0])
    no_cos = np.column_stack([(llh[:, 1] - lon0) * k, (llh[:, 0] - lat0) * k, llh[:, 2] - h0])
    err_flat = np.linalg.norm(flat - exact, axis=1)
    err_crude = np.linalg.norm(crude - exact, axis=1)
    err_nocos = np.linalg.norm(no_cos - exact, axis=1)
    fig, ax = figure(width=6.2, height=3.6)
    ax.loglog(d, err_nocos, color=SERIES[3], label="111 km per degree, no cos(latitude)")
    ax.loglog(d, err_crude, color=SERIES[1], label="111 km per degree, with cos(latitude)")
    ax.loglog(d, err_flat, color=SERIES[0], label=r"$N(\varphi_0)$, $M(\varphi_0)$ flat approximation (2.7)")
    ax.set_xlabel("distance from the origin [m]")
    ax.set_ylabel("error vs exact ENU [m]")
    ax.set_title("Local approximations against the exact ENU conversion")
    ax.legend()
    save(fig, "02_local_approximations.png")
    for dist in (100.0, 1000.0, 10000.0):
        i = np.argmin(abs(d - dist))
        print(f"  at {d[i]:8.0f} m: flat {err_flat[i]:.4f} m, 111km+cos {err_crude[i]:.3f} m, "
              f"no cos {err_nocos[i]:.1f} m")

    # 2.3 UTM point scale and grid convergence across zone 34 at the latitude of the data
    lon = np.linspace(18.0, 24.0, 121)
    u = geo.geodetic_to_utm(np.full_like(lon, lat0), lon, 34)
    fig, axes = figure(1, 2, width=8.0, height=3.0)
    axes[0].plot(u["easting"] / 1000, (u["scale"] - 1) * 1e6, color=SERIES[0])
    axes[0].axhline(0, color="0.5", lw=0.8)
    axes[0].set_xlabel("easting [km]")
    axes[0].set_ylabel("(k - 1) [ppm]")
    axes[0].set_title(f"UTM 34 point scale at {lat0:.2f} N")
    axes[1].plot(lon, u["convergence_deg"], color=SERIES[1])
    axes[1].set_xlabel("longitude [deg]")
    axes[1].set_ylabel(r"grid convergence $\gamma$ [deg]")
    axes[1].set_title("angle from true north to grid north")
    save(fig, "02_utm_scale.png")
    here = geo.geodetic_to_utm(lat0, lon0)
    print(f"  zone {here['zone']}: k = {here['scale']:.7f} ({(here['scale'] - 1) * 1e6:.0f} ppm), "
          f"gamma = {here['convergence_deg']:.4f} deg; k at the zone edge {u['scale'][0]:.6f}, "
          f"at the central meridian {u['scale'][60]:.6f}")


def static_and_walk():
    """ENU tracks of the two NMEA logs (fixes only), origin at the first fix, heights h = H + N_g."""
    tracks = {}
    for name in ["nmea_log1.nmea", "nmea_log2.nmea"]:
        epochs = nmea.parse_log_file(str(sf.data.path("gnss/" + name))).epochs
        a = sf.tables.epochs_to_arrays(epochs)
        v = a["fix_quality"] > 0
        h = a["alt_msl"] + a["geoid_sep"]
        ltp = geo.LocalTangentPlane(a["lat"][v][0], a["lon"][v][0], h[v][0])
        a["enu"] = np.full((len(v), 3), np.nan)
        a["enu"][v] = ltp.geodetic_to_enu(a["lat"][v], a["lon"][v], h[v])
        a["valid"] = v
        tracks[name] = (a, epochs)
    return tracks


def chapter3() -> None:
    print("chapter 3")
    tracks = static_and_walk()

    # 3.1 sky plot of the walk at one epoch
    a, epochs = tracks["nmea_log2.nmea"]
    sky = sf.tables.sky_per_epoch(epochs)
    k = 300
    s = sky[k]
    fig = plt_polar(s, f"nmea_log2, epoch {k}: {s['used'].sum()} used of {len(s['prn'])} in view")
    save(fig, "03_skyplot.png")
    d = gnss.compute_dop(s["azimuth_deg"][s["used"]], s["elevation_deg"][s["used"]])
    print(f"  epoch {k}: computed PDOP {d.pdop:.2f} HDOP {d.hdop:.2f} VDOP {d.vdop:.2f}; "
          f"reported {a['pdop'][k]:.2f} {a['hdop'][k]:.2f} {a['vdop'][k]:.2f}")

    for name in ["nmea_log1.nmea", "nmea_log2.nmea"]:
        an, ep = tracks[name]
        sk = sf.tables.sky_per_epoch(ep)
        close = total = 0
        for i, e in enumerate(ep):
            u = sk[i]["used"]
            if e.pdop is None or u.sum() < 4 or u.sum() != len(e.used_prns):
                continue
            dd = gnss.compute_dop(sk[i]["azimuth_deg"][u], sk[i]["elevation_deg"][u])
            total += 1
            close += abs(dd.pdop - e.pdop) < 0.1
        print(f"  {name}: reported PDOP within 0.1 of the GPS-only geometry in {close} of {total} epochs")

    # 3.2 DOP of a ring of n satellites plus one at the zenith, eq. (3.5)
    el = np.linspace(0, 75, 151)
    fig, ax = figure(width=6.0, height=3.4)
    for n, color in [(4, SERIES[0]), (8, SERIES[1])]:
        s_, c_ = np.sin(np.radians(el)), np.cos(np.radians(el))
        ax.plot(el, np.sqrt(4 / (n * c_**2)), color=color, label=f"HDOP, n = {n}")
        ax.plot(el, np.sqrt((n + 1) / (n * (1 - s_) ** 2)), color=color, ls="--", label=f"VDOP, n = {n}")
    ax.set_ylim(0, 8)
    ax.set_xlabel(r"elevation of the ring $\epsilon$ [deg]")
    ax.set_ylabel("DOP")
    ax.set_title("Ring of n satellites plus one at the zenith (3.5)")
    ax.legend(ncol=2)
    save(fig, "03_dop_ring.png")

    # 3.3 static log: scatter, 95 % ellipse, CEP circle
    a1, _ = tracks["nmea_log1.nmea"]
    v = a1["valid"] & (a1["fix_type"] == 3)
    en = a1["enu"][v][:, :2]
    st = gnss.precision_stats(en)
    ell = gnss.error_ellipse(st.covariance, 0.95)
    fig, ax = figure(width=5.0, height=4.6)
    ax.plot(en[:, 0] - st.mean[0], en[:, 1] - st.mean[1], ".-", ms=3, lw=0.4, color=SERIES[0], label="3D fixes")
    tt = np.linspace(0, 2 * np.pi, 200)
    xy = np.column_stack([ell.semi_major * np.cos(tt), ell.semi_minor * np.sin(tt)])
    R = np.array([[np.cos(ell.angle), -np.sin(ell.angle)], [np.sin(ell.angle), np.cos(ell.angle)]])
    xy = xy @ R.T
    ax.plot(xy[:, 0], xy[:, 1], color=SERIES[1], label="95 % ellipse (Gaussian)")
    ax.plot(st.cep50 * np.cos(tt), st.cep50 * np.sin(tt), color=SERIES[2], label="CEP (50 %, empirical)")
    ax.set_aspect("equal")
    ax.set_xlabel("east - mean [m]")
    ax.set_ylabel("north - mean [m]")
    ax.set_title("Standing still: nmea_log1, 3D fixes")
    ax.legend(loc="lower left")
    save(fig, "03_static_scatter.png")
    print(f"  static: n {st.count}, sigma_e {st.sigma_east:.2f} m, sigma_n {st.sigma_north:.2f} m, "
          f"DRMS {st.drms:.2f}, CEP {st.cep50:.2f}, R95 {st.r95:.2f}, 95% ellipse {ell.semi_major:.2f} x "
          f"{ell.semi_minor:.2f} m")

    # 3.4 displacement noise vs time lag on the static log (eq. 3.9)
    t = a1["t"][v]
    lags = np.array([1, 2, 3, 5, 7, 10, 15, 20, 30, 45, 60])
    sig = []
    for lag in lags:
        j = np.searchsorted(t, t + lag)
        ok = (j < len(t))
        ok[ok] &= np.abs(t[j[ok]] - t[ok] - lag) < 0.5
        dd = en[j[ok]] - en[ok]
        sig.append(np.sqrt(np.mean(dd**2)))  # mean over both axes = per-axis variance
    sig = np.array(sig)
    fig, ax = figure(width=6.0, height=3.4)
    ax.plot(lags, sig, "o-", color=SERIES[0], label="measured (static log)")
    ax.axhline(np.sqrt(2) * np.sqrt((st.sigma_east**2 + st.sigma_north**2) / 2), color=SERIES[1], ls="--",
               label=r"white-noise model $\sqrt{2}\,\sigma$")
    ax.set_xlabel(r"time lag $\tau$ [s]")
    ax.set_ylabel(r"$\sigma_d(\tau)$ per axis [m]")
    ax.set_title("Noise of a displacement over a time lag")
    ax.legend()
    save(fig, "03_displacement_noise.png")
    print("  sigma_d(tau): " + ", ".join(f"{l}s {x:.2f}" for l, x in zip(lags, sig)))

    # motion classification with the calibrated sigma_d(5 s)
    p = gnss.MotionParams()
    p.window_s = 5.0
    p.sigma_displacement_m = float(sig[lags == 5][0])
    for name in ["nmea_log1.nmea", "nmea_log2.nmea"]:
        an, _ = tracks[name]
        vv = an["valid"] & (an["fix_type"] == 3)
        for label, speed in [("position only", None), ("position + speed", an["speed"][vv])]:
            r = gnss.classify_motion(an["t"][vv], an["enu"][vv][:, :2], speed, p)
            print(f"  {name} {label}: moving in {100 * np.mean(r.moving):.0f} % of 3D epochs")


def plt_polar(s, title):
    import matplotlib.pyplot as plt
    from sensor_fusion.plotting import use_style

    use_style()
    fig = plt.figure(figsize=(4.6, 4.6), constrained_layout=True)
    ax = fig.add_subplot(projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    r = 90 - s["elevation_deg"]
    th = np.radians(s["azimuth_deg"])
    ax.scatter(th[s["used"]], r[s["used"]], s=60, color=SERIES[0], label="used in the fix", zorder=3)
    ax.scatter(th[~s["used"]], r[~s["used"]], s=60, facecolors="none", edgecolors=SERIES[1], label="in view, not used",
               zorder=3)
    for p, a_, b_ in zip(s["prn"], th, r):
        ax.annotate(str(p), (a_, b_), xytext=(6, 4), textcoords="offset points", fontsize=8, color="0.3")
    ax.set_rlim(0, 90)
    ax.set_rticks([30, 60, 90])
    ax.set_yticklabels(["el 60°", "30°", "0°"], color="0.45", fontsize=8)
    ax.set_rlabel_position(157)
    ax.set_title(title)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    return fig


CHAPTERS = {1: chapter1, 2: chapter2, 3: chapter3}

if __name__ == "__main__":
    wanted = [int(a) for a in sys.argv[1:]] or sorted(CHAPTERS)
    for c in wanted:
        CHAPTERS[c]()
