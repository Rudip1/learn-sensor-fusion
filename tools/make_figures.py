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


def demo_scan(noise: float = 0.02, seed: int = 1):
    """The 32-beam scan of chapter 4: sensor 1.8 m above the ground in the demo parking lot."""
    from sensor_fusion import sim

    T = np.eye(4)
    T[:3, 3] = [0.0, 0.0, 1.8]
    model = sim.LidarModel.spinning(32, -25.0, 10.0, 1024)
    model.range_noise = noise
    return sim.simulate_lidar(sim.demo_scene(), T, model, seed=seed)


def chapter4() -> None:
    from sensor_fusion import cloud

    print("chapter 4")
    scan = demo_scan()
    P, labels, true_n = scan["points"], scan["labels"], scan["normals"]
    print(f"  simulated scan: {len(P)} points, ground {np.mean(labels == -1):.2f}")
    fig, ax = figure(1, 3, width=13.5, height=4.0)
    # (a) the scan from above, ground grey, objects coloured by id
    g = labels == -1
    ax[0].scatter(P[g, 0], P[g, 1], s=0.2, color="0.75")
    ax[0].scatter(P[~g, 0], P[~g, 1], s=0.3, c=labels[~g] % 8, cmap="tab10")
    ax[0].set_aspect("equal")
    ax[0].set_xlim(-35, 35)
    ax[0].set_ylim(-25, 25)
    ax[0].set_title(f"simulated 32-beam scan ({len(P)} points)")
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("y [m]")
    # (b) voxel grid: points kept against the leaf size
    leaves = np.array([0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0])
    kept = [len(cloud.voxel_downsample(P, leaf)[0]) for leaf in leaves]
    for leaf, k in zip(leaves, kept):
        if leaf in (0.1, 0.5, 1.0):
            print(f"  voxel leaf {leaf} m: {k} points ({k / len(P):.2f})")
    ax[1].loglog(leaves, kept, "o-")
    ax[1].axhline(len(P), color="0.5", ls="--", lw=1)
    ax[1].set_xlabel("voxel edge [m]")
    ax[1].set_ylabel("points kept")
    ax[1].set_title("voxel-grid downsampling")
    # (c) normal error against the number of neighbours, 2 cm range noise
    for i, k in enumerate((5, 10, 20, 40)):
        n = cloud.estimate_normals(P, k, np.zeros(3))[0]
        ang = np.degrees(np.arccos(np.clip(np.abs(np.sum(n * true_n, axis=1)), 0, 1)))
        print(f"  normals k = {k}: median error {np.median(ang):.2f} deg, > 10 deg: {np.mean(ang > 10):.2f}")
        x = np.sort(ang)
        ax[2].plot(x, np.arange(1, len(x) + 1) / len(x), color=SERIES[i], label=f"k = {k}")
    ax[2].set_xscale("log")
    ax[2].set_xlim(0.05, 90)
    ax[2].set_xlabel("angle to the true normal [deg]")
    ax[2].set_ylabel("fraction of points")
    ax[2].set_title("PCA normals, 2 cm range noise")
    ax[2].legend()
    save(fig, "04_point_clouds.png")


def pose3(x, y, yaw, z):
    c, s = np.cos(yaw), np.sin(yaw)
    T = np.eye(4)
    T[:3, :3] = [[c, -s, 0], [s, c, 0], [0, 0, 1]]
    T[:3, 3] = [x, y, z]
    return T


def h2(x, y, a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, x], [s, c, y], [0, 0, 1.0]])


def planar_scans(poses2d, noise=0.01):
    """Planar scans (n, 2) of the demo scene from 2-D poses (x, y, yaw), scanner 1 m above the ground."""
    from sensor_fusion import sim

    model = sim.LidarModel.planar(720, 30.0)
    model.range_noise = noise
    scene = sim.demo_scene()
    return [sim.simulate_lidar(scene, pose3(x, y, a, 1.0), model, i)["points"][:, :2] for i, (x, y, a) in enumerate(poses2d)]


def recorded_scans():
    from sensor_fusion import cloud

    d = cloud.read_pcd(str(sf.data.path("parking_lot/run1_scans.pcd")))
    n = int(d["scan"].max()) + 1
    return [np.c_[d["x"][d["scan"] == k], d["y"][d["scan"] == k]] for k in range(n)]


def gnss_in_scan_window():
    """ENU positions (east, north) and times of the run-1 fixes during the recorded scans."""
    f = sf.io.load_navsatfix("parking_lot/run1_fix.csv")
    st = np.loadtxt(sf.data.path("parking_lot/run1_scans.csv"), delimiter=",", skiprows=1)
    ts, tf = st[:, 1] * 1e-9, f["stamp_ns"] * 1e-9
    ltp = geo.LocalTangentPlane(f["lat"][0], f["lon"][0], f["alt"][0])
    enu = ltp.geodetic_to_enu(f["lat"], f["lon"], f["alt"])
    w = (tf >= ts[0]) & (tf <= ts[-1])
    return ts, tf[w], enu[w, :2]


def chapter5() -> None:
    from sensor_fusion import reg, sim

    print("chapter 5")
    fig, ax = figure(1, 3, width=13.5, height=4.0)
    # (a) 3-D scan pair: convergence of the two metrics
    model = sim.LidarModel.spinning(32, -25.0, 10.0, 1024)
    model.range_noise = 0.01
    scene = sim.demo_scene()
    TA, TB = pose3(0, 0, 0, 1.8), pose3(1.0, 0.3, 0.05, 1.8)
    a, b = sim.simulate_lidar(scene, TA, model, 1), sim.simulate_lidar(scene, TB, model, 2)
    truth = np.linalg.inv(TA) @ TB
    for i, metric in enumerate(("point_to_point", "point_to_plane")):
        r = reg.icp(b["points"], a["points"], None, metric, 100, 2.0)
        E = np.linalg.inv(truth) @ r["T"]
        print(f"  3-D {metric}: {r['iterations']} iterations, error {np.linalg.norm(E[:3, 3]):.3f} m")
        ax[0].semilogy(np.arange(1, len(r["history"]) + 1), r["history"], color=SERIES[i],
                       label=f"{metric.replace('_', '-')}: {100 * np.linalg.norm(E[:3, 3]):.1f} cm off")
        objects = reg.icp(b["points"][b["labels"] >= 0], a["points"][a["labels"] >= 0], None, metric, 100, 2.0)
        E = np.linalg.inv(truth) @ objects["T"]
        print(f"  3-D {metric} without ground: error {np.linalg.norm(E[:3, 3]):.3f} m")
    ax[0].set_xlabel("iteration")
    ax[0].set_ylabel("RMS residual [m]")
    ax[0].set_title("ICP on two 32-beam scans (1.04 m apart)")
    ax[0].legend()
    # basin of convergence, 2-D
    A, B = planar_scans([(0, 0, 0), (0.5, 0.1, 0.03)])
    T2 = h2(0.5, 0.1, 0.03)
    rng = np.random.default_rng(0)
    for metric in ("point_to_point", "point_to_plane"):
        rates = []
        for off, yaw in ((1.0, 10), (2.0, 20), (3.0, 30)):
            ok = 0
            for _ in range(50):
                d = rng.normal(size=2)
                d *= off / np.linalg.norm(d)
                E = np.linalg.inv(T2) @ reg.icp(B, A, T2 @ h2(d[0], d[1], np.radians(yaw) * rng.choice([-1, 1])), metric, 100, 2.0)["T"]
                ok += np.linalg.norm(E[:2, 2]) < 0.05 and abs(np.arctan2(E[1, 0], E[0, 0])) < np.radians(0.5)
            rates.append(ok / 50)
        print(f"  basin {metric}: success from 1 m/10, 2 m/20, 3 m/30 deg: {rates}")
    # (b) simulated drive
    xs = np.linspace(-20, 20, 101)
    ys = 1.5 * np.sin(xs / 8)
    yaw = np.arctan(1.5 / 8 * np.cos(xs / 8))
    scans = planar_scans(zip(xs, ys, yaw))
    length = np.sum(np.hypot(np.diff(xs), np.diff(ys)))
    T0 = np.linalg.inv(h2(xs[0], ys[0], yaw[0]))
    rel = np.array([(T0 @ h2(x, y, a))[:2, 2] for x, y, a in zip(xs, ys, yaw)])
    ax[1].plot(rel[:, 0], rel[:, 1], "k", lw=2.5, label="truth")
    for i, metric in enumerate(("point_to_point", "point_to_plane")):
        r = reg.lidar_odometry(scans, metric, 50, 1.0)
        p = np.array([T[:2, 2] for T in r["poses"]])
        e = np.linalg.norm(p[-1] - rel[-1])
        print(f"  drive {length:.1f} m, {metric}: end error {e:.3f} m ({100 * e / length:.2f} %)")
        ax[1].plot(p[:, 0], p[:, 1], color=SERIES[i], label=f"{metric.replace('_', '-')}: {100 * e / length:.2f} % drift")
    ax[1].set_ylim(-4, 12)
    ax[1].set_xlabel("x [m]")
    ax[1].set_ylabel("y [m]")
    ax[1].set_title("odometry on a simulated 40 m drive")
    ax[1].legend(loc="upper left")
    # (c) recorded scans against GNSS
    real = recorded_scans()
    ts, tf, en = gnss_in_scan_window()
    ax[2].plot(en[:, 0] - en[0, 0], en[:, 1] - en[0, 1], ".", ms=3, color="0.5", label="GNSS fixes")
    for i, gate in enumerate((0.3, 0.5, 1.0)):
        r = reg.lidar_odometry(real, "point_to_plane", 50, gate)
        p = np.array([T[:2, 2] for T in r["poses"]])
        # odometry at the fix times, aligned rigidly to the fixes (the scanner's heading in ENU is unknown here)
        q = np.c_[np.interp(tf, ts, p[:, 0]), np.interp(tf, ts, p[:, 1])]
        A_ = reg.best_rigid(q, en)
        al = q @ A_[:2, :2].T + A_[:2, 2]
        rms = np.sqrt(np.mean(np.sum((al - en) ** 2, axis=1)))
        print(f"  recorded, gate {gate} m: net {np.linalg.norm(p[-1]):.1f} m, RMS to GNSS after alignment {rms:.2f} m")
        ax[2].plot(al[:, 0] - en[0, 0], al[:, 1] - en[0, 1], color=SERIES[i], label=f"gate {gate} m: RMS {rms:.1f} m")
    print(f"  GNSS: {len(en)} fixes, net {np.linalg.norm(en[-1] - en[0]):.1f} m")
    ax[2].set_xlim(-20, 20)
    ax[2].set_ylim(-36, 4)
    ax[2].set_xlabel("east [m]")
    ax[2].set_ylabel("north [m]")
    ax[2].set_title("recorded scans (61 s) against GNSS")
    ax[2].legend(fontsize=8, loc="upper left")
    save(fig, "05_registration.png")


def recorded_fusion_inputs():
    """Odometry (chapter 5, point-to-plane, gate 0.5 m) and GNSS of the recorded 61 s, and the initial state."""
    from sensor_fusion import reg

    real = recorded_scans()
    ts, tf, en = gnss_in_scan_window()
    f = sf.io.load_navsatfix("parking_lot/run1_fix.csv")
    w = np.isin(f["stamp_ns"] * 1e-9, tf)
    var = np.c_[f["var_e"][w], f["var_n"][w]]
    poses = reg.lidar_odometry(real, "point_to_plane", 50, 0.5)["poses"]
    odo = np.array([[T[0, 2], T[1, 2], np.arctan2(T[1, 0], T[0, 0])] for T in poses])
    k5 = np.argmax(np.linalg.norm(odo[:, :2], axis=1) > 5)
    g5 = np.argmax(np.linalg.norm(en - en[0], axis=1) > 5)
    yaw0 = np.arctan2(*(en[g5] - en[0])[::-1]) - np.arctan2(odo[k5, 1], odo[k5, 0])
    x0 = np.array([en[0, 0], en[0, 1], np.arctan2(np.sin(yaw0), np.cos(yaw0))])
    P0 = np.diag([var[0, 0], var[0, 1], np.radians(10) ** 2])
    return ts, odo, tf, en, var, x0, P0


def chapter6() -> None:
    from sensor_fusion import fusion, scenarios as sc

    print("chapter 6")
    noise = fusion.OdometryNoise(0.02, 0.05, 0.002)
    t, truth = sc.loop_drive(3.0)
    length = np.sum(np.linalg.norm(np.diff(truth[:, :2], axis=0), axis=1))
    P0 = np.diag([1.0, 1.0, np.radians(5) ** 2])
    print(f"  loop: {length:.0f} m, {t[-1]:.0f} s, odometry 10 Hz, GNSS 1 Hz sigma 1.5 m")
    cases = [("white GNSS", 0.0, {}), ("+ 5 % outliers", 0.0, {"outlier_fraction": 0.05}),
             ("+ 1 % odometry scale error", 0.01, {}), ("correlated GNSS (60 s)", 0.0, {"correlation_time": 60.0})]
    for name, scale, kw in cases:
        for gate in (0.99, 0.0):
            if gate == 0.0 and "outliers" not in name:
                continue
            st = []
            for seed in range(10):
                odo = sc.drifting_odometry(truth, noise, scale, rng=seed)
                tg, zg, _ = sc.gnss_fixes(t, truth, 1.0, 1.5, rng=100 + seed, **kw)
                tr = np.c_[np.interp(tg, t, truth[:, 0]), np.interp(tg, t, truth[:, 1])]
                r = fusion.fuse(list(t), odo, list(tg), zg, np.full((len(tg), 2), 2.25), truth[0], P0, noise, gate)
                e = r["x"][:, :2] - truth[:, :2]
                acc = np.array(r["accepted"])
                st.append([np.sqrt(np.mean(np.sum((zg - tr) ** 2, axis=1))), np.linalg.norm(odo[-1, :2] - truth[-1, :2]),
                           np.sqrt(np.mean(np.sum(e**2, axis=1))), np.linalg.norm(e, axis=1).max(),
                           np.mean(np.array(r["nis"])[acc]), np.sqrt(np.mean((e / r["sigma"][:, :2]) ** 2))])
            m = np.mean(st, axis=0)
            print(f"  {name:28s} gate {'on ' if gate else 'off'}: GNSS RMS {m[0]:.2f} m, odometry end {m[1]:.1f} m, "
                  f"fused RMS {m[2]:.2f} m (max {m[3]:.2f}), mean NIS {m[4]:.2f}, normalized error RMS {m[5]:.2f}")
    # one run of the outlier case for the figure
    odo = sc.drifting_odometry(truth, noise, 0.0, rng=0)
    tg, zg, out = sc.gnss_fixes(t, truth, 1.0, 1.5, outlier_fraction=0.05, rng=100)
    r = fusion.fuse(list(t), odo, list(tg), zg, np.full((len(tg), 2), 2.25), truth[0], P0, noise, 0.99)
    fig, ax = figure(1, 3, width=13.5, height=4.2)
    ax[0].plot(zg[~out, 0], zg[~out, 1], ".", ms=2, color="0.65", label="GNSS")
    ax[0].plot(zg[out, 0], zg[out, 1], "x", ms=5, color=SERIES[7], label="GNSS outliers")
    ax[0].plot(odo[:, 0], odo[:, 1], color=SERIES[1], lw=1, label="odometry")
    ax[0].plot(r["x"][:, 0], r["x"][:, 1], color=SERIES[0], lw=1.2, label="fused")
    ax[0].plot(truth[:, 0], truth[:, 1], "k--", lw=0.8, label="truth")
    ax[0].set_aspect("equal")
    ax[0].set_title("3 laps, 5 % GNSS outliers")
    ax[0].legend(fontsize=7, loc="center")
    e = np.linalg.norm(r["x"][:, :2] - truth[:, :2], axis=1)
    bound = 2 * np.sqrt(r["sigma"][:, 0] ** 2 + r["sigma"][:, 1] ** 2)
    tr = np.c_[np.interp(tg, t, truth[:, 0]), np.interp(tg, t, truth[:, 1])]
    ax[1].plot(tg, np.linalg.norm(zg - tr, axis=1), ".", ms=2, color="0.65", label="GNSS")
    ax[1].plot(t, np.linalg.norm(odo[:, :2] - truth[:, :2], axis=1), color=SERIES[1], label="odometry")
    ax[1].plot(t, e, color=SERIES[0], label="fused")
    ax[1].fill_between(t, 0, bound, color=SERIES[0], alpha=0.15, label="fused 2-sigma")
    ax[1].set_ylim(0, 8)
    ax[1].set_xlabel("time [s]")
    ax[1].set_ylabel("position error [m]")
    ax[1].legend(fontsize=7, loc="upper left")
    ax[1].set_title("errors against the truth")
    # recorded run
    ts, odo_r, tf, en, var, x0, P0r = recorded_fusion_inputs()
    for scale, color in ((1.0, SERIES[2]), (0.03, SERIES[0])):
        rr = fusion.fuse(list(ts), odo_r, list(tf), en, var, x0, P0r, noise, 0.99, scale)
        nis = np.mean(rr["nis"])
        print(f"  recorded: GNSS covariance x {scale}: mean NIS {nis:.2f}, accepted {np.mean(rr['accepted']):.2f}, "
              f"final heading sigma {np.degrees(rr['sigma'][-1, 2]):.2f} deg, heading {np.degrees(rr['x'][-1, 2]):.1f} deg")
        ax[2].plot(rr["x"][:, 0] - en[0, 0], rr["x"][:, 1] - en[0, 1], color=color, label=f"fused, R x {scale} (NIS {nis:.2f})")
    print(f"  recorded: reported GNSS sigma {np.sqrt(var).min():.2f}-{np.sqrt(var).max():.2f} m, initial heading {np.degrees(x0[2]):.1f} deg")
    ax[2].plot(en[:, 0] - en[0, 0], en[:, 1] - en[0, 1], ".", ms=2, color="0.6", label="GNSS fixes")
    ax[2].set_xlim(-20, 20)
    ax[2].set_ylim(-36, 4)
    ax[2].set_xlabel("east [m]")
    ax[2].set_ylabel("north [m]")
    ax[2].set_title("recorded run: odometry + GNSS")
    ax[2].legend(fontsize=7, loc="upper left")
    save(fig, "06_fusion.png")


def mapping_drive():
    """The chapter 7 drive: 16-beam LiDAR scans every 0.5 s along the lane, the true poses, and the scans' labels."""
    from sensor_fusion import mapping, scenarios as sc, sim

    t, truth = sc.lane_drive()
    T_vl = np.eye(4)
    T_vl[:3, 3] = [0, 0, 1.8]
    model = sim.LidarModel.spinning(16, -15.0, 15.0, 720, 60.0)
    model.range_noise = 0.02
    idx = np.arange(0, len(t), 5)
    scans = [sim.simulate_lidar(sim.demo_scene(), mapping.sensor_pose(truth[k], T_vl), model, int(k)) for k in idx]
    return t, truth, T_vl, idx, scans


def label_accuracy(points, labels, images, poses, cam, colors):
    """Fraction of coloured points whose own label equals the label of the pixel that coloured them."""
    from sensor_fusion import mapping

    cam_idx = np.array(colors["camera"])
    ok = n = 0
    for k, img in enumerate(images):
        sel = np.flatnonzero(cam_idx == k)
        if len(sel) == 0:
            continue
        uv = mapping.project(points[sel], poses[k], cam)["pixels"].astype(int)
        ok += np.sum(img["label"][uv[:, 1], uv[:, 0]] == labels[sel])
        n += len(sel)
    return n / len(points), ok / max(n, 1)


def chapter7() -> None:
    from sensor_fusion import fusion, mapping, scenarios as sc, sim

    print("chapter 7")
    t, truth, T_vl, idx, scans = mapping_drive()
    pts = [s["points"] for s in scans]
    print(f"  drive {np.sum(np.linalg.norm(np.diff(truth[:, :2], axis=0), axis=1)):.0f} m, {len(scans)} scans, {sum(len(p) for p in pts)} points")
    true_map = mapping.aggregate(pts, [mapping.sensor_pose(truth[k], T_vl) for k in idx])["points"]
    noise = fusion.OdometryNoise()
    odo = sc.drifting_odometry(truth, noise, 0.01, rng=1)
    tg, zg, _ = sc.gnss_fixes(t, truth, 1.0, 1.5, rng=2)
    fused = fusion.fuse(list(t), odo, list(tg), zg, np.full((len(tg), 2), 2.25), truth[0],
                        np.diag([1, 1, np.radians(5) ** 2]), noise)["x"]
    gx, gy = np.interp(t, tg, zg[:, 0]), np.interp(t, tg, zg[:, 1])
    gnss_only = np.c_[gx, gy, np.unwrap(np.arctan2(np.gradient(gy), np.gradient(gx)))]
    fig, ax = figure(1, 3, width=13.5, height=4.2)
    for i, (name, pose) in enumerate((("odometry", odo), ("fused", fused), ("GNSS only", gnss_only))):
        M = mapping.aggregate(pts, [mapping.sensor_pose(pose[k], T_vl) for k in idx])["points"]
        d = np.linalg.norm(M - true_map, axis=1)
        pe = np.linalg.norm(pose[idx, :2] - truth[idx, :2], axis=1)
        he = np.degrees(np.abs([fusion.wrap_angle(a) for a in pose[idx, 2] - truth[idx, 2]]))
        print(f"  {name:9s}: pose error RMS {np.sqrt(np.mean(pe**2)):.2f} m, heading {np.sqrt(np.mean(he**2)):.2f} deg; "
              f"map points displaced by {np.sqrt(np.mean(d**2)):.2f} m RMS, 95 % below {np.percentile(d, 95):.2f} m")
        x = np.sort(d)
        ax[2].plot(x, np.linspace(0, 1, len(x)), color=SERIES[i], label=name)
    ax[2].set_xscale("log")
    ax[2].set_xlim(0.01, 50)
    ax[2].set_xlabel("displacement of a map point [m]")
    ax[2].set_ylabel("fraction of points")
    ax[2].set_title("map error from pose error")
    ax[2].legend()
    # colouring with the true poses
    labels = np.concatenate([s["labels"] for s in scans])
    sub = np.random.default_rng(0).choice(len(true_map), 60000, replace=False)
    P, L = true_map[sub], labels[sub]
    cam = sim.Camera(np.array([[300.0, 0, 200], [0, 300, 150], [0, 0, 1]]), 400, 300)
    cams = np.arange(0, len(t), 25)
    Tc = [mapping.sensor_pose(truth[k], sc.forward_camera_mount()) for k in cams]
    images = [sim.render_camera(sim.demo_scene(), T, cam) for T in Tc]
    for occ in (False, True):
        c = mapping.colorize(P, [im["rgb"] for im in images], Tc, cam, occ)
        frac, acc = label_accuracy(P, L, images, Tc, cam, c)
        print(f"  colouring, occlusion test {'on ' if occ else 'off'}: {frac:.0%} of points coloured, {acc:.1%} from the right surface")
    for err in (0.5, 1.0, 2.0):
        Tbad = [mapping.sensor_pose(truth[k], sc.forward_camera_mount(yaw_error_deg=err)) for k in cams]
        imb = [sim.render_camera(sim.demo_scene(), T, cam) for T in Tbad]
        cb = mapping.colorize(P, [im["rgb"] for im in imb], Tc, cam, True)
        print(f"  camera yaw miscalibrated by {err} deg: {label_accuracy(P, L, imb, Tc, cam, cb)[1]:.1%} from the right surface")
    col = np.array(c["rgb"])
    seen = np.array(c["camera"]) >= 0
    ax[0].scatter(P[~seen, 0], P[~seen, 1], s=0.2, color="0.85")
    ax[0].scatter(P[seen, 0], P[seen, 1], s=0.4, c=np.clip(col[seen], 0, 1))
    ax[0].plot(truth[:, 0], truth[:, 1], "k", lw=1)
    ax[0].set_aspect("equal")
    ax[0].set_xlim(-35, 35)
    ax[0].set_ylim(-22, 22)
    ax[0].set_title(f"map coloured from {len(cams)} camera images")
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("y [m]")
    # one image with one scan projected
    k = 4
    img = images[k]
    sel = np.flatnonzero(np.isin(np.arange(len(true_map)), np.flatnonzero(np.repeat(np.arange(len(scans)), [len(p) for p in pts]) == cams[k] // 5)))
    pr = mapping.project(true_map[sel], Tc[k], cam)
    ins = np.array(pr["inside"])
    ax[1].imshow(np.clip(img["rgb"], 0, 1))
    ax[1].scatter(pr["pixels"][ins, 0], pr["pixels"][ins, 1], s=0.3, c=np.minimum(pr["depth"][ins], 40), cmap="turbo")
    ax[1].set_axis_off()
    ax[1].set_title("one scan projected into the camera (colour = depth)")
    save(fig, "07_map_building.png")


def detect_dynamic(drive, eps=0.5, min_points=3, voxel=0.5, ratio=0.3):
    """Chapter 8 pipeline on a scenarios.dynamic_drive(): RANSAC ground, DBSCAN of the rest, persistence filter.
    Returns per scan: world points, ground mask, cluster labels, dynamic mask, true mover mask."""
    from sensor_fusion import objects, scenarios as sc

    world, ground, labels = [], [], []
    for k, s in enumerate(drive):
        p, T = s["points"], s["T_world_sensor"]
        g = np.array(objects.segment_ground(p, 0.15, 500, 15.0, k)["ground"])
        lab = np.full(len(p), -1)
        lab[~g] = objects.dbscan(p[~g], eps, min_points)
        world.append(p @ T[:3, :3].T + T[:3, 3])
        ground.append(g)
        labels.append(lab)
    dyn = objects.dynamic_points(world, [l.tolist() for l in labels], voxel, ratio)
    movers = [np.isin(s["labels"], sc.MOVER_IDS) for s in drive]
    return world, ground, labels, [np.array(d) for d in dyn], movers


def chapter8() -> None:
    from sensor_fusion import objects, scenarios as sc

    print("chapter 8")
    for pitch in (0.0, 3.0):
        s = sc.dynamic_drive(0.0, 0.1, pitch_deg=pitch)[0]
        truth = s["labels"] == -1
        r = np.array(objects.segment_ground(s["points"])["ground"])
        h = np.array(objects.ground_by_height(s["points"], -1.8 + 0.15))
        print(f"  pitch {pitch} deg: RANSAC ground correct {np.mean(r == truth):.3f}, height threshold {np.mean(h == truth):.3f} "
              f"({np.sum(truth & ~h)} ground points missed)")
    fig, ax = figure(1, 3, width=13.5, height=4.2)
    results = {}
    for speed in (0.0, 2.0):
        drive = sc.dynamic_drive(speed)
        world, ground, labels, dyn, movers = detect_dynamic(drive)
        D, M = np.concatenate(dyn), np.concatenate(movers)
        tp = np.sum(D & M)
        G, L = np.concatenate(ground), np.concatenate(labels)
        print(f"  ego {speed} m/s: precision {tp / D.sum():.2f}, recall {tp / M.sum():.2f} ({M.sum()} mover points; "
              f"missed: {np.sum(M & G)} as ground, {np.sum(M & ~G & (L < 0))} as DBSCAN noise, {np.sum(M & ~D & (L >= 0))} in static clusters)")
        results[speed] = (world, ground, labels, dyn, movers)
    world, ground, labels, dyn, movers = results[0.0]
    # (a) one scan: clusters and the dynamic ones
    k = 15
    W, g, lab, d = world[k], ground[k], labels[k], dyn[k]
    ax[0].scatter(W[g, 0], W[g, 1], s=0.2, color="0.85")
    st = (lab >= 0) & ~d
    ax[0].scatter(W[st, 0], W[st, 1], s=0.5, c=lab[st] % 8, cmap="tab10")
    ax[0].scatter(W[d, 0], W[d, 1], s=3, color="k", label="dynamic")
    ax[0].scatter(W[(lab < 0) & ~g, 0], W[(lab < 0) & ~g, 1], s=0.5, color=SERIES[7], label="DBSCAN noise")
    ax[0].plot(-10, 0, "k^", ms=8)
    ax[0].set_aspect("equal")
    ax[0].set_xlim(-32, 32)
    ax[0].set_ylim(-20, 20)
    ax[0].set_title(f"scan at t = {1.5:.1f} s: clusters, dynamic in black")
    ax[0].legend(fontsize=7, loc="lower right")
    # (b) aggregated 3 s: ghost trails and what the filter removes
    allW = np.concatenate(world)
    allG, allD, allM = np.concatenate(ground), np.concatenate(dyn), np.concatenate(movers)
    keep = ~allG
    ax[1].scatter(allW[keep & ~allM, 0], allW[keep & ~allM, 1], s=0.2, color="0.6", label="static")
    ax[1].scatter(allW[allM & allD, 0], allW[allM & allD, 1], s=0.5, color=SERIES[0], label="moving, removed")
    ax[1].scatter(allW[allM & ~allD, 0], allW[allM & ~allD, 1], s=0.5, color=SERIES[7], label="moving, kept")
    ax[1].set_aspect("equal")
    ax[1].set_xlim(-32, 32)
    ax[1].set_ylim(-20, 20)
    ax[1].set_title("3 s of scans: the trails of the movers")
    ax[1].legend(fontsize=7, loc="lower right")
    # (c) persistence of static and moving points
    pers = np.concatenate(objects.persistence(world, 0.5))
    bins = np.linspace(0, 1, 31)
    ax[2].hist(pers[keep & ~allM], bins=bins, density=True, alpha=0.6, label="static objects")
    ax[2].hist(pers[allM], bins=bins, density=True, alpha=0.6, label="car and pedestrian")
    ax[2].axvline(0.3, color="k", ls="--", lw=1)
    ax[2].set_xlabel("persistence (fraction of the 30 scans)")
    ax[2].set_ylabel("density")
    ax[2].set_title("persistence, ego vehicle standing")
    ax[2].legend()
    print(f"  persistence: static median {np.median(pers[keep & ~allM]):.2f}, movers median {np.median(pers[allM]):.2f}")
    save(fig, "08_dynamic_objects.png")


CHAPTERS = {1: chapter1, 2: chapter2, 3: chapter3, 4: chapter4, 5: chapter5, 6: chapter6, 7: chapter7, 8: chapter8}

if __name__ == "__main__":
    wanted = [int(a) for a in sys.argv[1:]] or sorted(CHAPTERS)
    for c in wanted:
        CHAPTERS[c]()
