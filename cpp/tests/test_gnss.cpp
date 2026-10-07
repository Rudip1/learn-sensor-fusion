// Single-point positioning: exact recovery from noiseless pseudoranges and the least-squares optimality
// conditions with noise (1_theory/01_gnss_fundamentals.md, eqs. (1.1)-(1.5)).
#include <Eigen/Dense>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <random>

#include "sensor_fusion/gnss.hpp"

using Catch::Approx;
namespace gnss = sensor_fusion::gnss;

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;

// A receiver on a spherical Earth at 47.47 N, 19.06 E.
Eigen::Vector3d receiver() {
    const double lat = 47.47 * kDeg, lon = 19.06 * kDeg, R = 6'371'000.0;
    return R * Eigen::Vector3d(std::cos(lat) * std::cos(lon), std::cos(lat) * std::sin(lon), std::sin(lat));
}

Eigen::Matrix3Xd constellation(int n) {
    Eigen::VectorXd az(n), el(n);
    for (int i = 0; i < n; ++i) {
        az(i) = 2.0 * kPi * i / n + 0.3;
        el(i) = (15.0 + 60.0 * ((i * 37) % n) / n) * kDeg;
    }
    return gnss::place_satellites(receiver(), az, el);
}
}  // namespace

TEST_CASE("place_satellites reproduces azimuth, elevation and orbit radius", "[gnss][ch1]") {
    const Eigen::Vector3d r = receiver();
    Eigen::VectorXd az(3), el(3);
    az << 0.0, kPi / 2, 1.2;
    el << 10 * kDeg, 45 * kDeg, 80 * kDeg;
    const Eigen::Matrix3Xd s = gnss::place_satellites(r, az, el, 26'560'000.0);
    const Eigen::Vector3d up = r.normalized();
    const Eigen::Vector3d east = Eigen::Vector3d::UnitZ().cross(up).normalized();
    const Eigen::Vector3d north = up.cross(east);
    for (int i = 0; i < 3; ++i) {
        CHECK(s.col(i).norm() == Approx(26'560'000.0));
        const Eigen::Vector3d d = (s.col(i) - r).normalized();
        CHECK(std::asin(d.dot(up)) == Approx(el(i)));
        double a = std::atan2(d.dot(east), d.dot(north));
        if (a < 0) a += 2 * kPi;
        CHECK(a == Approx(az(i)).margin(1e-12));
    }
}

TEST_CASE("geometry matrix rows are [-u^T, 1] with unit u", "[gnss][ch1]") {
    const Eigen::Matrix3Xd s = constellation(6);
    const Eigen::MatrixX4d G = gnss::geometry_matrix(receiver(), s);
    for (int i = 0; i < 6; ++i) {
        CHECK(G.row(i).head<3>().norm() == Approx(1.0));
        CHECK(G(i, 3) == 1.0);
        // -u points from the satellite to the receiver: down-ish, so its dot with "up" is negative
        CHECK(G.row(i).head<3>().dot(receiver().normalized()) < 0.0);
    }
}

TEST_CASE("noiseless pseudoranges are inverted exactly", "[gnss][ch1]") {
    const double bias = 0.0123 * gnss::kSpeedOfLight;  // 12.3 ms receiver clock offset
    for (int n : {4, 5, 8, 12}) {
        const Eigen::Matrix3Xd s = constellation(n);
        const Eigen::VectorXd rho = gnss::geometric_ranges(receiver(), s).array() + bias;
        const auto sol = gnss::solve_position(s, rho);
        INFO("n = " << n);
        REQUIRE(sol.converged);
        CHECK((sol.position - receiver()).norm() < 1e-6);
        CHECK(sol.clock_bias_m == Approx(bias).margin(1e-6));
        CHECK(sol.residuals.cwiseAbs().maxCoeff() < 1e-6);
        CHECK(sol.iterations < 12);  // Gauss-Newton from the Earth's centre
    }
}

TEST_CASE("noisy solution satisfies the normal equations G^T W r = 0", "[gnss][ch1]") {
    const Eigen::Matrix3Xd s = constellation(9);
    std::mt19937 rng(7);
    std::normal_distribution<double> noise(0.0, 3.0);
    Eigen::VectorXd rho = gnss::geometric_ranges(receiver(), s).array() + 1000.0;
    for (int i = 0; i < rho.size(); ++i) rho(i) += noise(rng);
    Eigen::VectorXd w = Eigen::VectorXd::LinSpaced(9, 0.5, 2.0);
    const auto sol = gnss::solve_position(s, rho, w);
    REQUIRE(sol.converged);
    const Eigen::MatrixX4d G = gnss::geometry_matrix(sol.position, s);
    const Eigen::Vector4d grad = G.transpose() * w.asDiagonal() * sol.residuals;
    CHECK(grad.norm() < 1e-6);
    // cofactor is the inverse of the normal matrix, eq. (1.4)
    const Eigen::Matrix4d N = G.transpose() * w.asDiagonal() * G;
    CHECK((sol.cofactor * N - Eigen::Matrix4d::Identity()).norm() < 1e-9);
    CHECK((sol.position - receiver()).norm() < 30.0);
}

TEST_CASE("fewer than four satellites cannot be solved", "[gnss][ch1]") {
    const Eigen::Matrix3Xd s = constellation(3);
    const Eigen::VectorXd rho = gnss::geometric_ranges(receiver(), s);
    const auto sol = gnss::solve_position(s, rho);
    CHECK_FALSE(sol.converged);
    CHECK(sol.iterations == 0);
}

TEST_CASE("all satellites in one direction give a singular geometry", "[gnss][ch1]") {
    Eigen::VectorXd az = Eigen::VectorXd::Constant(5, 1.0);
    Eigen::VectorXd el = Eigen::VectorXd::Constant(5, 40 * kDeg);
    // identical lines of sight: G has rank 2 -> the solver must refuse
    const Eigen::Matrix3Xd s = gnss::place_satellites(receiver(), az, el);
    const Eigen::VectorXd rho = gnss::geometric_ranges(receiver(), s);
    gnss::SolverOptions opt;
    opt.initial_position = receiver() + Eigen::Vector3d(100, -50, 20);
    CHECK_FALSE(gnss::solve_position(s, rho, Eigen::VectorXd(), opt).converged);
}
