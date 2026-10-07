#include "sensor_fusion/gnss.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

namespace sensor_fusion::gnss {

Eigen::MatrixX4d geometry_matrix(const Eigen::Vector3d& receiver,
                                 const Eigen::Ref<const Eigen::Matrix3Xd>& satellites) {
    const Eigen::Index n = satellites.cols();
    Eigen::MatrixX4d G(n, 4);
    for (Eigen::Index i = 0; i < n; ++i) {
        const Eigen::Vector3d u = (satellites.col(i) - receiver).normalized();  // line of sight
        G.row(i) << -u.transpose(), 1.0;                                        // eq. (1.2)
    }
    return G;
}

Eigen::VectorXd geometric_ranges(const Eigen::Vector3d& receiver,
                                 const Eigen::Ref<const Eigen::Matrix3Xd>& satellites) {
    return (satellites.colwise() - receiver).colwise().norm().transpose();
}

PositionSolution solve_position(const Eigen::Ref<const Eigen::Matrix3Xd>& satellites,
                                const Eigen::Ref<const Eigen::VectorXd>& pseudoranges,
                                const Eigen::Ref<const Eigen::VectorXd>& weights, const SolverOptions& options) {
    const Eigen::Index n = satellites.cols();
    if (pseudoranges.size() != n) throw std::invalid_argument("solve_position: one pseudorange per satellite");
    if (weights.size() != 0 && weights.size() != n)
        throw std::invalid_argument("solve_position: weights must be empty or one per satellite");

    PositionSolution sol;
    sol.position = options.initial_position;
    sol.clock_bias_m = options.initial_clock_bias_m;
    if (n < 4) return sol;  // four unknowns need four equations

    const Eigen::VectorXd w = weights.size() == 0 ? Eigen::VectorXd::Ones(n) : Eigen::VectorXd(weights);
    for (int it = 0; it < options.max_iterations; ++it) {
        // Linearise about the current estimate, eq. (1.2).
        const Eigen::MatrixX4d G = geometry_matrix(sol.position, satellites);
        const Eigen::VectorXd predicted = geometric_ranges(sol.position, satellites).array() + sol.clock_bias_m;
        const Eigen::VectorXd dy = pseudoranges - predicted;

        // Normal equations, eq. (1.3).
        const Eigen::Matrix4d N = G.transpose() * w.asDiagonal() * G;
        const Eigen::LDLT<Eigen::Matrix4d> ldlt(N);
        if (ldlt.info() != Eigen::Success || ldlt.rcond() < 1e-12) return sol;  // singular geometry
        const Eigen::Vector4d dx = ldlt.solve(G.transpose() * w.asDiagonal() * dy);

        sol.position += dx.head<3>();
        sol.clock_bias_m += dx(3);
        sol.iterations = it + 1;
        if (dx.norm() < options.tolerance_m) {
            sol.converged = true;
            break;
        }
    }

    const Eigen::MatrixX4d G = geometry_matrix(sol.position, satellites);
    sol.cofactor = (G.transpose() * w.asDiagonal() * G).inverse();  // eq. (1.4)
    sol.residuals = pseudoranges - (geometric_ranges(sol.position, satellites).array() + sol.clock_bias_m).matrix();
    return sol;
}

Eigen::Matrix3Xd place_satellites(const Eigen::Vector3d& receiver, const Eigen::Ref<const Eigen::VectorXd>& azimuth,
                                  const Eigen::Ref<const Eigen::VectorXd>& elevation, double orbit_radius) {
    if (azimuth.size() != elevation.size()) throw std::invalid_argument("place_satellites: size mismatch");
    if (orbit_radius <= receiver.norm()) throw std::invalid_argument("place_satellites: orbit below receiver");

    // Local east-north-up basis with "up" along the geocentric radius.
    const Eigen::Vector3d up = receiver.normalized();
    Eigen::Vector3d east = Eigen::Vector3d::UnitZ().cross(up);
    if (east.norm() < 1e-9) east = Eigen::Vector3d::UnitY();  // at a pole any east will do
    east.normalize();
    const Eigen::Vector3d north = up.cross(east);

    Eigen::Matrix3Xd sats(3, azimuth.size());
    for (Eigen::Index i = 0; i < azimuth.size(); ++i) {
        const double ce = std::cos(elevation(i));
        const Eigen::Vector3d d =
            ce * std::sin(azimuth(i)) * east + ce * std::cos(azimuth(i)) * north + std::sin(elevation(i)) * up;
        // |r + t d| = R  ->  t^2 + 2 t (r.d) + |r|^2 - R^2 = 0, positive root. Eq. (1.5).
        const double rd = receiver.dot(d);
        const double t = -rd + std::sqrt(rd * rd - receiver.squaredNorm() + orbit_radius * orbit_radius);
        sats.col(i) = receiver + t * d;
    }
    return sats;
}

}  // namespace sensor_fusion::gnss
