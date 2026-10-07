#pragma once
/// @file gnss.hpp
/// Single-point positioning from pseudoranges (Gauss-Newton least squares) and a satellite-geometry helper.
/// Theory: 1_theory/01_gnss_fundamentals.md, section "How a fix is computed".

#include <Eigen/Core>

namespace sensor_fusion::gnss {

/// Speed of light in vacuum [m/s] (exact, SI definition).
constexpr double kSpeedOfLight = 299792458.0;

struct SolverOptions {
    int max_iterations = 20;
    double tolerance_m = 1e-4;                                      ///< stop when |dx| falls below this
    Eigen::Vector3d initial_position = Eigen::Vector3d::Zero();    ///< ECEF start; Earth's centre works
    double initial_clock_bias_m = 0.0;
};

struct PositionSolution {
    Eigen::Vector3d position = Eigen::Vector3d::Zero();  ///< receiver position, ECEF [m]
    double clock_bias_m = 0.0;                           ///< b = c * dt_r [m]
    Eigen::Matrix4d cofactor = Eigen::Matrix4d::Zero();  ///< (G^T W G)^-1, eq. (1.4)
    Eigen::VectorXd residuals;                           ///< rho - rho_hat at the solution [m]
    int iterations = 0;
    bool converged = false;
};

/// Rows g_i^T = [-u_i^T, 1] of the geometry matrix, u_i the unit line of sight receiver -> satellite i.
/// Eq. (1.2). @p satellites is 3 x n (ECEF).
Eigen::MatrixX4d geometry_matrix(const Eigen::Vector3d& receiver, const Eigen::Ref<const Eigen::Matrix3Xd>& satellites);

/// Geometric ranges |s_i - r| [m].
Eigen::VectorXd geometric_ranges(const Eigen::Vector3d& receiver,
                                 const Eigen::Ref<const Eigen::Matrix3Xd>& satellites);

/// Solve rho_i = |s_i - r| + b for (r, b) by Gauss-Newton, eqs. (1.1)-(1.3).
/// @param weights per-measurement weights w_i = 1 / sigma_i^2; pass an empty vector for unit weights.
/// Needs at least four satellites; with fewer, or a singular geometry, `converged` is false.
PositionSolution solve_position(const Eigen::Ref<const Eigen::Matrix3Xd>& satellites,
                                const Eigen::Ref<const Eigen::VectorXd>& pseudoranges,
                                const Eigen::Ref<const Eigen::VectorXd>& weights = Eigen::VectorXd(),
                                const SolverOptions& options = SolverOptions());

/// Place satellites at azimuth/elevation (radians, as seen from @p receiver) on a sphere of radius
/// @p orbit_radius around the Earth's centre. Uses the local vertical r/|r| (spherical Earth), which is
/// enough to build test geometries. Eq. (1.5).
Eigen::Matrix3Xd place_satellites(const Eigen::Vector3d& receiver, const Eigen::Ref<const Eigen::VectorXd>& azimuth,
                                  const Eigen::Ref<const Eigen::VectorXd>& elevation,
                                  double orbit_radius = 26'560'000.0);

}  // namespace sensor_fusion::gnss
