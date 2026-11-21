"""Synthetic drives with exact ground truth for the fusion chapters: a true trajectory, odometry that drifts
(eq. 6.3 noise plus a scale error) and GNSS fixes with white or time-correlated noise and outliers."""

from __future__ import annotations

import numpy as np

from ._core import fusion


def loop_drive(laps: float = 1.0, radius: float = 40.0, speed: float = 5.0, odom_rate: float = 10.0):
    """True poses (x, y, yaw) of a vehicle driving a rounded-rectangle loop (two 2*radius straights joined by
    half circles) at constant speed. Returns (times, poses (n, 3))."""
    straight, arc = 2 * radius, np.pi * radius
    lap = 2 * straight + 2 * arc
    n = int(laps * lap / speed * odom_rate) + 1
    s = np.arange(n) * speed / odom_rate
    poses = np.zeros((n, 3))
    for i, si in enumerate(s % lap):
        if si < straight:                                   # bottom straight, eastwards
            poses[i] = [si - radius, -radius, 0.0]
        elif si < straight + arc:                           # right half circle
            a = (si - straight) / radius
            poses[i] = [radius + radius * np.sin(a), -radius * np.cos(a), a]
        elif si < 2 * straight + arc:                       # top straight, westwards
            poses[i] = [radius - (si - straight - arc), radius, np.pi]
        else:                                               # left half circle
            a = (si - 2 * straight - arc) / radius
            poses[i] = [-radius - radius * np.sin(a), radius * np.cos(a), np.pi + a]
    poses[:, 2] = np.vectorize(fusion.wrap_angle)(poses[:, 2])
    return s / speed, poses


def drifting_odometry(truth: np.ndarray, noise: "fusion.OdometryNoise", scale_error: float = 0.01, rng=None):
    """Odometry poses from the true increments, each scaled by (1 + scale_error) and perturbed with the
    covariance noise.covariance(delta) (eq. 6.3). Starts at the true first pose."""
    rng = np.random.default_rng(rng)
    out = [truth[0].copy()]
    for k in range(1, len(truth)):
        d = np.array(fusion.between(truth[k - 1], truth[k]))
        d[:2] *= 1 + scale_error
        d = rng.multivariate_normal(d, noise.covariance(d))
        out.append(np.array(fusion.compose(out[-1], d)))
    return np.array(out)


def gnss_fixes(times: np.ndarray, truth: np.ndarray, rate: float = 1.0, sigma: float = 1.5,
               correlation_time: float = 0.0, outlier_fraction: float = 0.0, outlier_size: float = 20.0, rng=None):
    """GNSS fixes at `rate` Hz sampled from the true positions: noise of standard deviation `sigma` per axis,
    white (correlation_time = 0) or first-order Gauss-Markov with that time constant, and a fraction of
    outliers displaced by `outlier_size` metres in a random direction. Returns (t, positions (m, 2), outlier mask)."""
    rng = np.random.default_rng(rng)
    t = np.arange(times[0], times[-1] + 1e-9, 1.0 / rate)
    xy = np.c_[np.interp(t, times, truth[:, 0]), np.interp(t, times, truth[:, 1])]
    e = np.zeros_like(xy)
    if correlation_time > 0:
        phi = np.exp(-1.0 / (rate * correlation_time))
        e[0] = rng.normal(0, sigma, 2)
        for i in range(1, len(t)):
            e[i] = phi * e[i - 1] + rng.normal(0, sigma * np.sqrt(1 - phi**2), 2)
    else:
        e = rng.normal(0, sigma, xy.shape)
    out = rng.random(len(t)) < outlier_fraction
    ang = rng.uniform(0, 2 * np.pi, out.sum())
    e[out] += outlier_size * np.c_[np.cos(ang), np.sin(ang)]
    return t, xy + e, out


def lane_drive(duration: float = 25.0, speed: float = 2.0, odom_rate: float = 10.0):
    """True poses (x, y, yaw) of a slow drive along the lane of sim.demo_scene(), from x = -25 m eastwards with a
    gentle 0.8 m sway. Returns (times, poses (n, 3))."""
    t = np.arange(0.0, duration + 1e-9, 1.0 / odom_rate)
    x = -25.0 + speed * t
    y = 0.8 * np.sin(x / 10)
    yaw = np.arctan(0.08 * np.cos(x / 10))
    return t, np.c_[x, y, yaw]


def forward_camera_mount(x: float = 0.5, height: float = 1.5, yaw_error_deg: float = 0.0) -> np.ndarray:
    """T_vehicle_camera of a camera looking along the vehicle's x axis (camera z forward, x right, y down),
    optionally rotated by a yaw (calibration) error."""
    a = np.radians(yaw_error_deg)
    Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    T = np.eye(4)
    T[:3, :3] = Rz @ np.column_stack([(0, -1, 0), (0, 0, -1), (1, 0, 0)])
    T[:3, 3] = [x, 0, height]
    return T
