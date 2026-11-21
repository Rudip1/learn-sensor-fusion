// Chapter 5 tests: the closed-form alignment, ICP with both metrics on simulated scans with known motion
// (including the bias of point-to-point ICP on the ground rings of a 3-D LiDAR), and LiDAR odometry on a
// simulated drive and on the recorded scans.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <random>

#include "sensor_fusion/registration.hpp"
#include "sensor_fusion/sim.hpp"

using namespace sensor_fusion;
using Catch::Approx;

namespace {

Eigen::Isometry3d pose3(double x, double y, double yaw, double z) {
    Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
    T.linear() = Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    T.translation() = Eigen::Vector3d(x, y, z);
    return T;
}

reg::Rigid<2> pose2(double x, double y, double yaw) {
    reg::Rigid<2> T = reg::Rigid<2>::Identity();
    T.linear() = Eigen::Rotation2Dd(yaw).toRotationMatrix();
    T.translation() = Eigen::Vector2d(x, y);
    return T;
}

cloud::Points<2> planar_scan(const sim::Scene& scene, const reg::Rigid<2>& T, std::uint32_t seed) {
    auto m = sim::LidarModel::planar(720, 30.0);
    m.range_noise = 0.01;
    const auto s = sim::simulate_lidar(scene, pose3(T.translation().x(), T.translation().y(),
                                                    Eigen::Rotation2Dd(T.linear()).angle(), 1.0), m, seed);
    return s.points.topRows<2>();
}

/// Translation [m] and rotation [rad] of truth^-1 * estimate.
template <int Dim>
std::pair<double, double> error(const reg::Rigid<Dim>& truth, const reg::Rigid<Dim>& est) {
    const reg::Rigid<Dim> E = truth.inverse() * est;
    double angle = 0.0;
    if constexpr (Dim == 2)
        angle = std::abs(Eigen::Rotation2Dd(E.linear()).angle());
    else
        angle = Eigen::AngleAxisd(E.linear()).angle();
    return {E.translation().norm(), angle};
}

}  // namespace

TEST_CASE("best_rigid recovers a rigid transform exactly and never a reflection", "[reg][ch5]") {
    std::mt19937 gen(3);
    std::normal_distribution<double> G(0, 1);
    cloud::Points<3> a(3, 20);
    for (int i = 0; i < 20; ++i) a.col(i) << G(gen), G(gen), G(gen);
    reg::Rigid<3> T = reg::Rigid<3>::Identity();
    T.linear() = Eigen::AngleAxisd(0.7, Eigen::Vector3d(1, 2, 3).normalized()).toRotationMatrix();
    T.translation() = Eigen::Vector3d(1, -2, 0.5);
    const auto e = error<3>(T, reg::best_rigid<3>(a, cloud::transform<3>(T, a)));
    CHECK(e.first < 1e-12);
    CHECK(e.second < 1e-7);
    // 2-D
    cloud::Points<2> b = a.topRows<2>();
    const reg::Rigid<2> T2 = pose2(3, 1, -2.0);
    CHECK(error<2>(T2, reg::best_rigid<2>(b, cloud::transform<2>(T2, b))).first < 1e-12);
    // coplanar points (z = 0) mirrored in z: the best proper rotation is returned, det = +1
    cloud::Points<3> flat = a;
    flat.row(2).setZero();
    cloud::Points<3> mirrored = flat;
    mirrored.row(0) *= -1;
    CHECK(reg::best_rigid<3>(flat, mirrored).linear().determinant() == Approx(1.0));
}

TEST_CASE("ICP on identical clouds stops at once", "[reg][ch5]") {
    const auto s = planar_scan(sim::demo_scene(), pose2(0, 0, 0), 1);
    const auto r = reg::icp<2>(s, s, reg::Rigid<2>::Identity());
    CHECK(r.converged);
    CHECK(r.iterations == 1);
    CHECK(r.rmse < 1e-12);
}

TEST_CASE("3-D ICP: point-to-plane is accurate; point-to-point is pulled by the ground rings", "[reg][ch5]") {
    const sim::Scene scene = sim::demo_scene();
    auto m = sim::LidarModel::spinning(32, -25.0 * M_PI / 180, 10.0 * M_PI / 180, 1024);
    m.range_noise = 0.01;
    const Eigen::Isometry3d TA = pose3(0, 0, 0, 1.8), TB = pose3(1.0, 0.3, 0.05, 1.8);
    const auto a = sim::simulate_lidar(scene, TA, m, 1), b = sim::simulate_lidar(scene, TB, m, 2);
    reg::Rigid<3> truth;
    truth.matrix() = (TA.inverse() * TB).matrix();
    reg::IcpOptions o;
    o.max_iterations = 100;
    o.max_correspondence = 2.0;
    o.metric = reg::Metric::PointToPlane;
    const auto plane = reg::icp<3>(b.points, a.points, reg::Rigid<3>::Identity(), o);
    const auto ep = error<3>(truth, plane.T);
    CHECK(ep.first < 0.03);
    CHECK(ep.second < 0.005);
    o.trim = 0.8;  // dropping the 20 % worst pairs (parts seen in one scan only) makes it exact to the noise
    CHECK(error<3>(truth, reg::icp<3>(b.points, a.points, reg::Rigid<3>::Identity(), o).T).first < 0.005);
    // point-to-point: the rings on the ground move with the sensor and hold the estimate back
    o.trim = 1.0;
    o.metric = reg::Metric::PointToPoint;
    const auto point = reg::icp<3>(b.points, a.points, reg::Rigid<3>::Identity(), o);
    CHECK(error<3>(truth, point.T).first > 0.1);
    CHECK(point.iterations > 3 * plane.iterations);
    // ...and is accurate once the ground is removed (chapter 8 removes it without labels)
    auto objects = [](const sim::LidarScan& s) {
        std::vector<int> keep;
        for (int i = 0; i < s.labels.size(); ++i)
            if (s.labels(i) >= 0) keep.push_back(i);
        cloud::Points<3> p(3, Eigen::Index(keep.size()));
        for (std::size_t k = 0; k < keep.size(); ++k) p.col(Eigen::Index(k)) = s.points.col(keep[k]);
        return p;
    };
    CHECK(error<3>(truth, reg::icp<3>(objects(b), objects(a), reg::Rigid<3>::Identity(), o).T).first < 0.03);
}

TEST_CASE("2-D ICP: both metrics converge from 1 m and 10 degrees; point-to-plane in fewer iterations", "[reg][ch5]") {
    const sim::Scene scene = sim::demo_scene();
    const reg::Rigid<2> truth = pose2(0.5, 0.1, 0.03);
    const auto A = planar_scan(scene, pose2(0, 0, 0), 1), B = planar_scan(scene, truth, 2);
    reg::IcpOptions o;
    o.max_iterations = 100;
    o.max_correspondence = 2.0;
    int it[2] = {0, 0};
    for (int m = 0; m < 2; ++m) {
        o.metric = m == 0 ? reg::Metric::PointToPoint : reg::Metric::PointToPlane;
        for (double a : {-1.0, 1.0}) {
            const auto r = reg::icp<2>(B, A, truth * pose2(0.7 * a, -0.7 * a, a * 10 * M_PI / 180), o);
            const auto e = error<2>(truth, r.T);
            CHECK(e.first < 0.05);
            CHECK(e.second < 0.5 * M_PI / 180);
            it[m] += r.iterations;
        }
    }
    CHECK(it[1] < it[0]);
}

TEST_CASE("LiDAR odometry on a simulated drive: drift below 1 % with point-to-plane", "[reg][ch5]") {
    const sim::Scene scene = sim::demo_scene();
    std::vector<reg::Rigid<2>> truth;
    std::vector<cloud::Points<2>> scans;
    const int n = 100;
    double length = 0;
    for (int k = 0; k <= n; ++k) {
        const double x = -20.0 + 40.0 * k / n, y = 1.5 * std::sin(x / 8);
        const double yaw = std::atan(1.5 / 8 * std::cos(x / 8));
        truth.push_back(pose2(x, y, yaw));
        if (k > 0) length += (truth[std::size_t(k)].translation() - truth[std::size_t(k - 1)].translation()).norm();
        scans.push_back(planar_scan(scene, truth.back(), std::uint32_t(k)));
    }
    reg::OdometryOptions o;
    o.icp.max_correspondence = 1.0;
    double end_error[2];
    for (int m = 0; m < 2; ++m) {
        o.icp.metric = m == 0 ? reg::Metric::PointToPoint : reg::Metric::PointToPlane;
        const auto r = reg::lidar_odometry<2>(scans, o);
        REQUIRE(r.poses.size() == truth.size());
        end_error[m] = error<2>(truth.front().inverse() * truth.back(), r.poses.back()).first;
    }
    CHECK(end_error[1] < 0.01 * length);
    CHECK(end_error[0] > 3 * end_error[1]);
}

TEST_CASE("LiDAR odometry on the recorded planar scans", "[reg][ch5][data]") {
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
    REQUIRE(scans.size() == 484);
    reg::OdometryOptions o;
    o.icp.metric = reg::Metric::PointToPlane;
    o.icp.max_correspondence = 0.5;
    const auto r = reg::lidar_odometry<2>(scans, o);
    // GNSS over the same 61 s: 33.7 m between the first and the last fix (chapter 6)
    const double net = r.poses.back().translation().norm();
    CHECK(net > 31.0);
    CHECK(net < 36.0);
    for (std::size_t k = 1; k < r.poses.size(); ++k)  // at 8 Hz no step is longer than 0.5 m
        CHECK((r.poses[k].translation() - r.poses[k - 1].translation()).norm() < 0.5);
}
