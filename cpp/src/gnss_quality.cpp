#include "sensor_fusion/gnss_quality.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "sensor_fusion/geodesy.hpp"

namespace sensor_fusion::gnss {

namespace {

Dop dop_from_enu_cofactor(const Eigen::Matrix4d& Q) {
    Dop d;
    d.cofactor_enu = Q;
    d.hdop = std::sqrt(Q(0, 0) + Q(1, 1));  // eq. (3.3)
    d.vdop = std::sqrt(Q(2, 2));
    d.pdop = std::sqrt(Q(0, 0) + Q(1, 1) + Q(2, 2));
    d.tdop = std::sqrt(Q(3, 3));
    d.gdop = std::sqrt(Q.trace());
    d.valid = true;
    return d;
}

double quantile_sorted(const std::vector<double>& v, double q) {
    if (v.empty()) return 0.0;
    const double pos = q * static_cast<double>(v.size() - 1);
    const auto i = static_cast<std::size_t>(std::floor(pos));
    const double frac = pos - static_cast<double>(i);
    return i + 1 < v.size() ? v[i] * (1 - frac) + v[i + 1] * frac : v[i];
}

}  // namespace

Eigen::MatrixX4d geometry_matrix_enu(const Eigen::Ref<const Eigen::VectorXd>& azimuth,
                                     const Eigen::Ref<const Eigen::VectorXd>& elevation) {
    if (azimuth.size() != elevation.size()) throw std::invalid_argument("azimuth and elevation sizes differ");
    Eigen::MatrixX4d G(azimuth.size(), 4);
    for (Eigen::Index i = 0; i < azimuth.size(); ++i) {
        const double ce = std::cos(elevation(i));
        // line of sight u = [cos(el) sin(az), cos(el) cos(az), sin(el)] in ENU; row = [-u^T, 1], eq. (3.1)
        G.row(i) << -ce * std::sin(azimuth(i)), -ce * std::cos(azimuth(i)), -std::sin(elevation(i)), 1.0;
    }
    return G;
}

Dop compute_dop(const Eigen::Ref<const Eigen::VectorXd>& azimuth, const Eigen::Ref<const Eigen::VectorXd>& elevation,
                const Eigen::Ref<const Eigen::VectorXd>& weights) {
    if (weights.size() != 0 && weights.size() != azimuth.size())
        throw std::invalid_argument("weights must be empty or one per satellite");
    if (azimuth.size() < 4) return Dop{};
    const Eigen::MatrixX4d G = geometry_matrix_enu(azimuth, elevation);
    const Eigen::VectorXd w = weights.size() == 0 ? Eigen::VectorXd::Ones(azimuth.size()) : Eigen::VectorXd(weights);
    const Eigen::Matrix4d N = G.transpose() * w.asDiagonal() * G;  // eq. (3.2)
    const Eigen::FullPivLU<Eigen::Matrix4d> lu(N);
    if (lu.rank() < 4 || lu.rcond() < 1e-12) return Dop{};
    return dop_from_enu_cofactor(lu.inverse());
}

Dop dop_from_ecef_cofactor(const Eigen::Matrix4d& cofactor_ecef, double lat, double lon) {
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topLeftCorner<3, 3>() = geo::rotation_ecef_enu(lat, lon).transpose();  // ECEF -> ENU
    return dop_from_enu_cofactor(T * cofactor_ecef * T.transpose());       // eq. (3.4)
}

Eigen::VectorXd elevation_weights(const Eigen::Ref<const Eigen::VectorXd>& elevation) {
    return elevation.array().sin().square().matrix();
}

ErrorEllipse error_ellipse(const Eigen::Matrix2d& cov, double probability) {
    if (!(probability > 0.0 && probability < 1.0)) throw std::invalid_argument("probability must be in (0, 1)");
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es(cov);
    const Eigen::Vector2d lambda = es.eigenvalues().cwiseMax(0.0);  // ascending
    // The squared Mahalanobis radius of a 2-D Gaussian is chi^2 with 2 dof: P(r^2 < c) = 1 - exp(-c/2).
    const double c = -2.0 * std::log(1.0 - probability);  // eq. (3.7)
    ErrorEllipse e;
    e.semi_major = std::sqrt(c * lambda(1));
    e.semi_minor = std::sqrt(c * lambda(0));
    const Eigen::Vector2d major = es.eigenvectors().col(1);
    e.angle = std::atan2(major.y(), major.x());
    return e;
}

PrecisionStats precision_stats(const Eigen::Ref<const Eigen::Matrix2Xd>& en) {
    PrecisionStats s;
    s.count = static_cast<int>(en.cols());
    if (s.count == 0) return s;
    s.mean = en.rowwise().mean();
    const Eigen::Matrix2Xd d = en.colwise() - s.mean;
    if (s.count > 1) s.covariance = d * d.transpose() / static_cast<double>(s.count - 1);
    s.sigma_east = std::sqrt(s.covariance(0, 0));
    s.sigma_north = std::sqrt(s.covariance(1, 1));
    s.drms = std::sqrt(s.covariance.trace());  // eq. (3.6)
    s.two_drms = 2.0 * s.drms;
    std::vector<double> r(static_cast<std::size_t>(s.count));
    for (int i = 0; i < s.count; ++i) r[static_cast<std::size_t>(i)] = d.col(i).norm();
    std::sort(r.begin(), r.end());
    s.cep50 = quantile_sorted(r, 0.5);
    s.r95 = quantile_sorted(r, 0.95);
    return s;
}

std::string gga_quality_name(int quality) {
    switch (quality) {
        case 0: return "invalid";
        case 1: return "GPS (standard positioning)";
        case 2: return "differential (DGPS / SBAS)";
        case 3: return "PPS";
        case 4: return "RTK fixed";
        case 5: return "RTK float";
        case 6: return "dead reckoning";
        case 7: return "manual input";
        case 8: return "simulator";
        default: return "unknown";
    }
}

std::string gsa_fix_type_name(int fix_type) {
    switch (fix_type) {
        case 1: return "no fix";
        case 2: return "2D";
        case 3: return "3D";
        default: return "unknown";
    }
}

std::string navsat_status_name(int status) {
    switch (status) {
        case -1: return "STATUS_NO_FIX";
        case 0: return "STATUS_FIX";
        case 1: return "STATUS_SBAS_FIX";
        case 2: return "STATUS_GBAS_FIX";
        default: return "unknown";
    }
}

MotionResult classify_motion(const Eigen::Ref<const Eigen::VectorXd>& t, const Eigen::Ref<const Eigen::Matrix2Xd>& en,
                             const Eigen::Ref<const Eigen::VectorXd>& speed, const MotionParams& params) {
    const Eigen::Index n = t.size();
    if (en.cols() != n) throw std::invalid_argument("classify_motion: one position per time stamp");
    if (speed.size() != 0 && speed.size() != n) throw std::invalid_argument("classify_motion: speed size");
    const double nan = std::numeric_limits<double>::quiet_NaN();

    MotionResult r;
    r.test_statistic = Eigen::VectorXd::Constant(n, nan);
    r.displacement_speed = Eigen::VectorXd::Constant(n, nan);
    std::vector<bool> raw(static_cast<std::size_t>(n), false);

    const double var_diff = params.sigma_displacement_m * params.sigma_displacement_m;
    Eigen::Index j = 0;  // earliest epoch at least window_s before i
    for (Eigen::Index i = 0; i < n; ++i) {
        while (j + 1 < i && t(j + 1) <= t(i) - params.window_s) ++j;
        bool moving = false;
        if (i > 0 && t(i) - t(j) >= params.window_s * 0.999) {
            const Eigen::Vector2d dp = en.col(i) - en.col(j);
            r.test_statistic(i) = dp.squaredNorm() / var_diff;  // eq. (3.8)
            r.displacement_speed(i) = dp.norm() / (t(i) - t(j));
            moving = r.test_statistic(i) > params.chi2_threshold;
        }
        if (speed.size() == n && std::isfinite(speed(i)) && speed(i) > params.speed_threshold_mps) moving = true;
        raw[static_cast<std::size_t>(i)] = moving;
    }

    // Hysteresis: a change of state is accepted only if it lasts min_run epochs.
    r.moving.assign(static_cast<std::size_t>(n), false);
    if (n == 0) return r;
    bool state = raw[0];
    for (Eigen::Index i = 0; i < n; ++i) {
        if (raw[static_cast<std::size_t>(i)] != state) {
            Eigen::Index k = i;
            while (k < n && raw[static_cast<std::size_t>(k)] != state && k - i < params.min_run) ++k;
            if (k - i >= params.min_run || k == n) state = !state;
        }
        r.moving[static_cast<std::size_t>(i)] = state;
    }
    return r;
}

}  // namespace sensor_fusion::gnss
