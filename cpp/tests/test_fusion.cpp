// Chapter 6 tests: pose composition, the EKF prediction (Jacobians against finite differences), the update
// (against the closed-form Kalman gain), the chi-square gate, consistency of the filter on a simulated loop
// with white GNSS noise, rejection of outliers, and the recorded run (pessimistic reported covariances).
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <limits>
#include <array>
#include <fstream>
#include <random>
#include <sstream>

#include "sensor_fusion/cloud.hpp"
#include "sensor_fusion/fusion.hpp"
#include "sensor_fusion/geodesy.hpp"
#include "sensor_fusion/registration.hpp"

using namespace sensor_fusion;
using namespace sensor_fusion::fusion;
using Catch::Approx;

TEST_CASE("compose and between are inverse; angles wrap", "[fusion][ch6]") {
    const Eigen::Vector3d a(1.0, -2.0, 3.0), d(0.4, 0.1, 0.5);
    const Eigen::Vector3d b = compose(a, d);
    CHECK((between(a, b) - d).norm() < 1e-12);
    CHECK(b(2) == Approx(wrap_angle(3.5)));
    CHECK(wrap_angle(3.5) == Approx(3.5 - 2 * M_PI));
    CHECK(wrap_angle(-M_PI) == Approx(M_PI));
}

TEST_CASE("EKF prediction propagates the covariance with the Jacobians of eq. (6.1)", "[fusion][ch6]") {
    const Eigen::Vector3d x0(2.0, 1.0, 0.7), d(0.8, -0.2, 0.1);
    const Eigen::Matrix3d P0 = Eigen::Vector3d(0.3, 0.2, 0.05).asDiagonal();
    const Eigen::Matrix3d Q = Eigen::Vector3d(0.01, 0.02, 0.003).asDiagonal();
    PoseEkf ekf(x0, P0);
    ekf.predict(d, Q);
    CHECK((ekf.state() - compose(x0, d)).norm() < 1e-12);
    // numerical Jacobians of compose with respect to the pose and to the increment
    Eigen::Matrix3d F, G;
    const double h = 1e-7;
    for (int i = 0; i < 3; ++i) {
        const Eigen::Vector3d e = h * Eigen::Vector3d::Unit(i);
        F.col(i) = (compose(x0 + e, d) - compose(x0 - e, d)) / (2 * h);
        G.col(i) = (compose(x0, d + e) - compose(x0, d - e)) / (2 * h);
    }
    const Eigen::Matrix3d P = F * P0 * F.transpose() + G * Q * G.transpose();
    CHECK((ekf.covariance() - P).norm() < 1e-7);
}

TEST_CASE("EKF update equals the Kalman gain in closed form; the gate rejects outliers", "[fusion][ch6]") {
    PoseEkf ekf(Eigen::Vector3d(0, 0, 0.2), Eigen::Vector3d(4.0, 4.0, 0.01).asDiagonal());
    const auto u = ekf.update_position(Eigen::Vector2d(1.0, -2.0), Eigen::Matrix2d::Identity(), chi2_2dof(0.99));
    CHECK(u.accepted);
    CHECK(u.nis == Approx((1.0 + 4.0) / 5.0));
    // gain 4 / (4 + 1) = 0.8 on each axis, uncorrelated heading untouched
    CHECK(ekf.state()(0) == Approx(0.8));
    CHECK(ekf.state()(1) == Approx(-1.6));
    CHECK(ekf.state()(2) == Approx(0.2));
    CHECK(ekf.covariance()(0, 0) == Approx(0.8));
    // a fix 20 m away: NIS = 400 / (0.8 + 1) >> 9.21, rejected, state unchanged
    const Eigen::Vector3d before = ekf.state();
    const auto o = ekf.update_position(Eigen::Vector2d(20.8, -1.6), Eigen::Matrix2d::Identity(), chi2_2dof(0.99));
    CHECK_FALSE(o.accepted);
    CHECK(o.nis == Approx(400.0 / 1.8));
    CHECK(ekf.state() == before);
    CHECK(chi2_2dof(0.99) == Approx(9.2103).epsilon(1e-4));
}

namespace {

struct Loop {
    std::vector<double> t;
    std::vector<Eigen::Vector3d> truth, odom;
    std::vector<bool> outlier;  ///< per GNSS fix
    FusionInput in;
};

/// Three laps of a 40 m-radius circle at 5 m/s, odometry at 10 Hz with the noise model itself (no unmodelled
/// error), GNSS at 1 Hz with white noise of 1.5 m and a fraction of 20 m outliers.
Loop make_loop(double outliers, unsigned seed, const OdometryNoise& noise) {
    std::mt19937 gen(seed);
    std::normal_distribution<double> N(0, 1);
    std::uniform_real_distribution<double> U(0, 1);
    Loop L;
    const double r = 40, v = 5, dt = 0.1;
    const int n = int(3 * 2 * M_PI * r / v / dt);
    for (int k = 0; k <= n; ++k) {
        const double a = v * k * dt / r;
        L.t.push_back(k * dt);
        L.truth.emplace_back(r * std::sin(a), r - r * std::cos(a), wrap_angle(a));
    }
    L.odom.push_back(L.truth[0]);
    for (int k = 1; k <= n; ++k) {
        Eigen::Vector3d d = between(L.truth[std::size_t(k - 1)], L.truth[std::size_t(k)]);
        const Eigen::Vector3d s = noise.covariance(d).diagonal().cwiseSqrt();
        d += Eigen::Vector3d(s(0) * N(gen), s(1) * N(gen), s(2) * N(gen));
        L.odom.push_back(compose(L.odom.back(), d));
    }
    L.in.odom_times = L.t;
    L.in.odom_poses = L.odom;
    for (int k = 0; k <= n; k += 10) {
        Eigen::Vector2d z = L.truth[std::size_t(k)].head<2>() + 1.5 * Eigen::Vector2d(N(gen), N(gen));
        L.outlier.push_back(U(gen) < outliers);
        if (L.outlier.back()) {
            const double phi = 2 * M_PI * U(gen);
            z += 20.0 * Eigen::Vector2d(std::cos(phi), std::sin(phi));
        }
        L.in.gnss_times.push_back(L.t[std::size_t(k)]);
        L.in.gnss_positions.push_back(z);
        L.in.gnss_covariances.push_back(2.25 * Eigen::Matrix2d::Identity());
    }
    return L;
}

}  // namespace

TEST_CASE("fusion on a simulated loop: better than both sensors and consistent", "[fusion][ch6]") {
    FusionOptions o;
    const Loop L = make_loop(0.0, 7, o.odometry);
    const auto r = fuse(L.in, L.truth[0], Eigen::Vector3d(1, 1, 0.01).asDiagonal(), o);
    REQUIRE(r.states.size() == L.truth.size());
    double se = 0, sz = 0, so = 0, sg = 0, nis = 0;
    for (std::size_t k = 0; k < L.truth.size(); ++k) {
        const Eigen::Vector2d e = r.states[k].head<2>() - L.truth[k].head<2>();
        se += e.squaredNorm();
        sz += e.dot(r.covariances[k].topLeftCorner<2, 2>().ldlt().solve(e));  // NEES of the position
        so += (L.odom[k].head<2>() - L.truth[k].head<2>()).squaredNorm();
    }
    for (std::size_t g = 0; g < L.in.gnss_times.size(); ++g) {
        sg += (L.in.gnss_positions[g] - L.truth[g * 10].head<2>()).squaredNorm();
        nis += r.nis[g];
    }
    const double n = double(L.truth.size()), m = double(L.in.gnss_times.size());
    const double fused = std::sqrt(se / n), gnss = std::sqrt(sg / m), odom = std::sqrt(so / n);
    CHECK(fused < 0.6 * gnss);
    CHECK(fused < odom);
    CHECK(sz / n > 2 * 0.7);  // mean NEES of a consistent 2-D estimate is 2
    CHECK(sz / n < 2 * 1.4);
    CHECK(nis / m > 2 * 0.8);
    CHECK(nis / m < 2 * 1.25);
}

TEST_CASE("the chi-square gate rejects GNSS outliers", "[fusion][ch6]") {
    FusionOptions o;
    const Loop L = make_loop(0.05, 8, o.odometry);
    auto max_error = [&](const FusionResult& r) {
        double e = 0;
        for (std::size_t k = 0; k < L.truth.size(); ++k) e = std::max(e, (r.states[k] - L.truth[k]).head<2>().norm());
        return e;
    };
    const auto gated = fuse(L.in, L.truth[0], Eigen::Vector3d(1, 1, 0.01).asDiagonal(), o);
    o.gate_probability = 0.0;
    const auto open = fuse(L.in, L.truth[0], Eigen::Vector3d(1, 1, 0.01).asDiagonal(), o);
    int outliers = 0, caught = 0, false_alarms = 0;
    for (std::size_t g = 0; g < gated.accepted.size(); ++g) {
        outliers += L.outlier[g];
        caught += L.outlier[g] && !gated.accepted[g];
        false_alarms += !L.outlier[g] && !gated.accepted[g];
    }
    REQUIRE(outliers > 0);
    CHECK(caught == outliers);  // 20 m against 1.5 m noise: every outlier fails the 99 % gate
    CHECK(false_alarms <= 0.03 * double(gated.accepted.size()));
    CHECK(max_error(gated) < 0.6 * max_error(open));
}

TEST_CASE("recorded run: the receiver's covariances are far too pessimistic", "[fusion][ch6][data]") {
    // odometry from the recorded scans (chapter 5) and the GNSS fixes of the same 61 s
    const auto d = cloud::read_pcd(std::string(SF_DATA_DIR) + "/parking_lot/run1_scans.pcd");
    const Eigen::VectorXd sx = d.field("x"), sy = d.field("y"), sid = d.field("scan");
    const int n = int(sid.maxCoeff()) + 1;
    std::vector<std::vector<Eigen::Vector2d>> pts(static_cast<std::size_t>(n));
    for (Eigen::Index i = 0; i < sx.size(); ++i) pts[std::size_t(sid(i))].emplace_back(sx(i), sy(i));
    std::vector<cloud::Points<2>> scans;
    for (const auto& p : pts) {
        cloud::Points<2> s(2, Eigen::Index(p.size()));
        for (std::size_t k = 0; k < p.size(); ++k) s.col(Eigen::Index(k)) = p[k];
        scans.push_back(s);
    }
    reg::OdometryOptions ro;
    ro.icp.metric = reg::Metric::PointToPlane;
    ro.icp.max_correspondence = 0.5;
    const auto odo = reg::lidar_odometry<2>(scans, ro);
    FusionInput in;
    std::ifstream stamps(std::string(SF_DATA_DIR) + "/parking_lot/run1_scans.csv");
    std::string line;
    std::getline(stamps, line);
    while (std::getline(stamps, line)) in.odom_times.push_back(std::stod(line.substr(line.find(',') + 1)) * 1e-9);
    REQUIRE(in.odom_times.size() == odo.poses.size());
    for (const auto& T : odo.poses)
        in.odom_poses.emplace_back(T.translation().x(), T.translation().y(), Eigen::Rotation2Dd(T.linear()).angle());
    std::ifstream fixes(std::string(SF_DATA_DIR) + "/parking_lot/run1_fix.csv");
    std::getline(fixes, line);
    std::vector<std::array<double, 7>> rows;
    while (std::getline(fixes, line)) {
        std::array<double, 7> v{};
        std::stringstream ss(line);
        std::string cell;
        // stamp_ns, lat, lon, alt, status, service, var_e, var_n, var_u
        std::vector<double> c;
        while (std::getline(ss, cell, ',')) c.push_back(std::stod(cell));
        v = {c[0] * 1e-9, c[1], c[2], c[3], c[6], c[7], 0};
        rows.push_back(v);
    }
    auto geod = [](const std::array<double, 7>& v) { return geo::Geodetic{geo::deg2rad(v[1]), geo::deg2rad(v[2]), v[3]}; };
    const geo::LocalTangentPlane ltp(geod(rows[0]));
    for (const auto& v : rows) {
        if (v[0] < in.odom_times.front() || v[0] > in.odom_times.back()) continue;
        in.gnss_times.push_back(v[0]);
        in.gnss_positions.push_back(ltp.geodetic_to_enu(geod(v)).head<2>());
        in.gnss_covariances.push_back(Eigen::Vector2d(v[4], v[5]).asDiagonal());
    }
    REQUIRE(in.gnss_times.size() == 243);
    // initial heading: GNSS direction over the first 5 m minus the odometry direction over its first 5 m
    std::size_t g5 = 0, k5 = 0;
    while ((in.gnss_positions[g5] - in.gnss_positions[0]).norm() < 5) ++g5;
    while (in.odom_poses[k5].head<2>().norm() < 5) ++k5;
    const Eigen::Vector2d dg = in.gnss_positions[g5] - in.gnss_positions[0];
    const double yaw0 = std::atan2(dg.y(), dg.x()) - std::atan2(in.odom_poses[k5].y(), in.odom_poses[k5].x());
    const Eigen::Vector3d x0(in.gnss_positions[0].x(), in.gnss_positions[0].y(), wrap_angle(yaw0));
    FusionOptions o;
    const auto r = fuse(in, x0, Eigen::Vector3d(in.gnss_covariances[0](0, 0), in.gnss_covariances[0](1, 1), 0.03).asDiagonal(), o);
    double nis = 0;
    for (double v : r.nis) nis += v;
    CHECK(nis / double(r.nis.size()) < 0.1);  // a consistent filter would average 2
}
