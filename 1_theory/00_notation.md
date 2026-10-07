# 0 · Notation

One notation is used through the whole module. Symbols are defined here once; a chapter repeats a
definition only where it is first needed.

## Scalars, vectors, matrices

| symbol | meaning |
|---|---|
| $a$, $\alpha$ | scalar (italic) |
| $\mathbf{p}$, $\mathbf{x}$ | column vector (bold lower case) |
| $\mathbf{R}$, $\mathbf{P}$ | matrix (bold upper case) |
| $\mathbf{I}_n$ | $n\times n$ identity |
| $\lVert\mathbf{v}\rVert$ | Euclidean norm |
| $\hat{x}$ | estimate of $x$ |
| $\check{x}$ | predicted (a priori) value of $x$, before a measurement update |
| $\oplus$ | bitwise exclusive or (chapter 1 only) |

## Frames and transforms

| symbol | frame |
|---|---|
| $\{E\}$ | ECEF, Earth-centred Earth-fixed (WGS-84) |
| $\{N\}$ | local ENU (east, north, up) tangent plane at a chosen origin; the "world" or "map" frame |
| $\{B\}$ | body (vehicle, ego) frame: $x$ forward, $y$ left, $z$ up |
| $\{L\}$ | LiDAR sensor frame |
| $\{C\}$ | camera frame: $z$ along the optical axis, $x$ right, $y$ down |

- ${}^{A}\mathbf{p}$ is a point expressed in frame $\{A\}$.
- $\mathbf{T}_{AB}\in SE(3)$ maps coordinates from $\{B\}$ to $\{A\}$: ${}^{A}\mathbf{p} = \mathbf{R}_{AB}\,{}^{B}\mathbf{p} + \mathbf{t}_{AB}$,
  i.e. $\mathbf{T}_{AB}$ is the pose of $\{B\}$ in $\{A\}$. Chains read right to left:
  $\mathbf{T}_{AC}=\mathbf{T}_{AB}\mathbf{T}_{BC}$.
- Homogeneous coordinates: $\tilde{\mathbf{p}} = [\mathbf{p}^\top\ 1]^\top$,
  $\mathbf{T} = \begin{bmatrix}\mathbf{R} & \mathbf{t}\\ \mathbf{0}^\top & 1\end{bmatrix}$.
- Planar poses: $(x, y, \theta)\in SE(2)$, $\theta$ the heading (yaw), counter-clockwise from the $x$ axis.

## GNSS (chapters 1–3)

| symbol | meaning |
|---|---|
| $\varphi,\ \lambda,\ h$ | geodetic latitude, longitude, ellipsoidal height |
| $H$, $N_g$ | orthometric (mean-sea-level) height, geoid undulation; $h = H + N_g$ |
| $a,\ f,\ b,\ e^2$ | ellipsoid semi-major axis, flattening, semi-minor axis, first eccentricity squared |
| $N(\varphi),\ M(\varphi)$ | prime-vertical and meridian radii of curvature |
| $\mathbf{r}$ | receiver position (ECEF) |
| $\mathbf{s}_i$ | position of satellite $i$ (ECEF) |
| $\rho_i$ | pseudorange to satellite $i$ |
| $c$ | speed of light, $299\,792\,458$ m/s |
| $\delta t_r$, $b = c\,\delta t_r$ | receiver clock offset, the same expressed in metres ("clock bias") |
| $\mathbf{u}_i$ | unit line-of-sight vector from receiver to satellite $i$ |
| $\mathbf{G}$ | geometry (design) matrix, rows $[-\mathbf{u}_i^\top\ \ 1]$ |
| $\mathbf{Q}$ | cofactor matrix $(\mathbf{G}^\top\mathbf{G})^{-1}$ |
| $\alpha,\ \epsilon$ | azimuth (clockwise from north) and elevation of a satellite |
| $\sigma_\rho$ | user-equivalent range error (pseudorange noise, 1σ) |

## Point clouds and registration (chapters 4–5, 7–8)

| symbol | meaning |
|---|---|
| $\mathcal{P}=\{\mathbf{p}_i\}_{i=1}^{n}$ | source point cloud |
| $\mathcal{Q}=\{\mathbf{q}_j\}_{j=1}^{m}$ | target (reference) point cloud |
| $\mathbf{n}_i$ | unit surface normal at $\mathbf{p}_i$ |
| $\ell$ | voxel edge length |
| $k$ | number of neighbours |
| $\bar{\mathbf{p}}$, $\boldsymbol{\Sigma}$ | centroid and covariance of a point set |
| $d_{\max}$ | maximum correspondence distance |
| $\varepsilon$, $m_{\min}$ | DBSCAN neighbourhood radius and minimum neighbour count |

## Estimation (chapter 6)

| symbol | meaning |
|---|---|
| $\mathbf{x}_k$, $\mathbf{P}_k$ | state and its covariance at step $k$ |
| $\mathbf{u}_k$ | control / odometry input |
| $\mathbf{F}_k$, $\mathbf{L}_k$ | Jacobians of the motion model w.r.t. state and input |
| $\mathbf{Q}_k$, $\mathbf{R}_k$ | process (input) and measurement noise covariances |
| $\mathbf{z}_k$, $\mathbf{H}_k$ | measurement and its Jacobian |
| $\boldsymbol{\nu}_k$, $\mathbf{S}_k$ | innovation and innovation covariance |
| $\mathbf{K}_k$ | Kalman gain |
| $d^2$ | squared Mahalanobis distance $\boldsymbol{\nu}^\top\mathbf{S}^{-1}\boldsymbol{\nu}$ |
| $\chi^2_{m}(p)$ | $p$-quantile of the chi-square distribution with $m$ degrees of freedom |

## Cameras (chapter 7)

| symbol | meaning |
|---|---|
| $\mathbf{K}_c$ | camera intrinsic matrix (subscript $c$ keeps it apart from the Kalman gain) |
| $f_x, f_y, c_x, c_y$ | focal lengths and principal point in pixels |
| $(u, v)$ | pixel coordinates, $u$ to the right, $v$ down |

## Code conventions

- C++ stores a point set as a $3\times n$ (or $2\times n$) Eigen matrix, one point per **column**.
  The Python module exposes the same data as $(n, 3)$ numpy arrays, one point per **row**, which is what
  numpy, SciPy and Open3D expect.
- Angles are radians in C++ unless a name ends in `_deg`. Distances are metres.
- A C++ comment `// eq. (3.4)` refers to equation (3.4) of the theory chapter with the same number.
