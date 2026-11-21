#include "sensor_fusion/fusion.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace sensor_fusion::fusion {

double wrap_angle(double a) {
    a = std::fmod(a + M_PI, 2 * M_PI);
    if (a <= 0) a += 2 * M_PI;
    return a - M_PI;
}

Eigen::Vector3d compose(const Eigen::Vector3d& p, const Eigen::Vector3d& d) {
    const double c = std::cos(p(2)), s = std::sin(p(2));
    return {p(0) + c * d(0) - s * d(1), p(1) + s * d(0) + c * d(1), wrap_angle(p(2) + d(2))};  // (6.1)
}

Eigen::Vector3d between(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    const double c = std::cos(a(2)), s = std::sin(a(2));
    const double dx = b(0) - a(0), dy = b(1) - a(1);
    return {c * dx + s * dy, -s * dx + c * dy, wrap_angle(b(2) - a(2))};
}

Eigen::Matrix3d OdometryNoise::covariance(const Eigen::Vector3d& d) const {
    const double dist = d.head<2>().norm();
    const double st = k_trans * dist + floor_trans;
    const double sr = k_rot * std::abs(d(2)) + k_rot_per_m * dist + floor_rot;  // (6.3)
    return Eigen::Vector3d(st * st, st * st, sr * sr).asDiagonal();
}

void PoseEkf::predict(const Eigen::Vector3d& d, const Eigen::Matrix3d& Q) {
    const double c = std::cos(x_(2)), s = std::sin(x_(2));
    Eigen::Matrix3d F = Eigen::Matrix3d::Identity();  // d f / d x
    F(0, 2) = -s * d(0) - c * d(1);
    F(1, 2) = c * d(0) - s * d(1);
    Eigen::Matrix3d G = Eigen::Matrix3d::Identity();  // d f / d u
    G.topLeftCorner<2, 2>() << c, -s, s, c;
    x_ = compose(x_, d);
    P_ = F * P_ * F.transpose() + G * Q * G.transpose();  // (6.2)
}

PoseEkf::Update PoseEkf::update_position(const Eigen::Vector2d& z, const Eigen::Matrix2d& R, double gate) {
    Update u;
    u.innovation = z - x_.head<2>();                               // (6.4), H = [I 0]
    const Eigen::Matrix2d S = P_.topLeftCorner<2, 2>() + R;
    u.nis = u.innovation.dot(S.ldlt().solve(u.innovation));        // (6.5)
    if (!(u.nis <= gate)) return u;
    const Eigen::Matrix<double, 3, 2> K = P_.leftCols<2>() * S.inverse();  // (6.6)
    x_ += K * u.innovation;
    x_(2) = wrap_angle(x_(2));
    // Joseph form: P = (I - K H) P (I - K H)^T + K R K^T, symmetric and positive semidefinite
    Eigen::Matrix<double, 2, 3> H = Eigen::Matrix<double, 2, 3>::Zero();
    H.leftCols<2>().setIdentity();
    const Eigen::Matrix3d IKH = Eigen::Matrix3d::Identity() - K * H;
    P_ = IKH * P_ * IKH.transpose() + K * R * K.transpose();
    u.accepted = true;
    return u;
}

double chi2_2dof(double p) { return -2.0 * std::log(1.0 - p); }

FusionResult fuse(const FusionInput& in, const Eigen::Vector3d& x0, const Eigen::Matrix3d& P0, const FusionOptions& o) {
    if (in.odom_times.size() != in.odom_poses.size() || in.odom_poses.empty())
        throw std::invalid_argument("fuse: odometry times and poses differ in length");
    if (in.gnss_times.size() != in.gnss_positions.size() || in.gnss_times.size() != in.gnss_covariances.size())
        throw std::invalid_argument("fuse: GNSS inputs differ in length");
    const double gate = (o.gate_probability > 0 && o.gate_probability < 1) ? chi2_2dof(o.gate_probability)
                                                                          : std::numeric_limits<double>::infinity();
    FusionResult r;
    PoseEkf ekf(x0, P0);
    std::size_t g = 0;
    auto absorb = [&](double t) {  // all fixes up to time t
        while (g < in.gnss_times.size() && in.gnss_times[g] <= t) {
            const auto u = ekf.update_position(in.gnss_positions[g], o.gnss_scale * in.gnss_covariances[g], gate);
            r.nis.push_back(u.nis);
            r.accepted.push_back(u.accepted);
            ++g;
        }
    };
    // fixes before the first odometry pose are skipped (no state yet)
    while (g < in.gnss_times.size() && in.gnss_times[g] < in.odom_times[0]) ++g;
    for (std::size_t k = 0; k < in.odom_poses.size(); ++k) {
        if (k > 0) {
            const Eigen::Vector3d d = between(in.odom_poses[k - 1], in.odom_poses[k]);
            ekf.predict(d, o.odometry.covariance(d));
        }
        const double t_next = k + 1 < in.odom_times.size() ? in.odom_times[k + 1] : std::numeric_limits<double>::infinity();
        // a fix is applied at the last odometry time before it (at most one odometry period late)
        absorb(std::nextafter(t_next, -std::numeric_limits<double>::infinity()));
        r.times.push_back(in.odom_times[k]);
        r.states.push_back(ekf.state());
        r.covariances.push_back(ekf.covariance());
    }
    return r;
}

}  // namespace sensor_fusion::fusion
