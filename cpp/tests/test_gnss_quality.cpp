// DOP against the closed form of a symmetric constellation (1_theory/03_gnss_quality.md, eq. 3.5), the
// ECEF route of chapter 1, error-ellipse probabilities and the motion classifier on synthetic tracks.
#include <Eigen/Dense>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <random>

#include "sensor_fusion/geodesy.hpp"
#include "sensor_fusion/gnss.hpp"
#include "sensor_fusion/gnss_quality.hpp"

using Catch::Approx;
namespace gnss = sensor_fusion::gnss;
namespace geo = sensor_fusion::geo;

namespace {
constexpr double kPi = 3.14159265358979323846;
}

TEST_CASE("DOP of a ring of n satellites plus one at the zenith: closed form (3.5)", "[dop][ch3]") {
    for (int n : {3, 4, 6, 10}) {
        for (double el_deg : {5.0, 20.0, 45.0}) {
            const double el = el_deg * kPi / 180.0, s = std::sin(el), c = std::cos(el);
            Eigen::VectorXd az(n + 1), elev(n + 1);
            for (int i = 0; i < n; ++i) {
                az(i) = 2 * kPi * i / n + 0.1;
                elev(i) = el;
            }
            az(n) = 0.0;
            elev(n) = kPi / 2;
            const gnss::Dop d = gnss::compute_dop(az, elev);
            INFO("n = " << n << ", elevation " << el_deg);
            REQUIRE(d.valid);
            CHECK(d.hdop == Approx(std::sqrt(4.0 / (n * c * c))));
            CHECK(d.vdop == Approx(std::sqrt((n + 1.0) / (n * (1 - s) * (1 - s)))));
            CHECK(d.tdop == Approx(std::sqrt((n * s * s + 1.0) / (n * (1 - s) * (1 - s)))));
            CHECK(d.pdop == Approx(std::sqrt(d.hdop * d.hdop + d.vdop * d.vdop)));
            CHECK(d.gdop == Approx(std::sqrt(d.pdop * d.pdop + d.tdop * d.tdop)));
        }
    }
}

TEST_CASE("DOP is invariant to a rotation of the sky about the vertical", "[dop][ch3]") {
    Eigen::VectorXd az(5), el(5);
    az << 0.3, 1.4, 2.6, 4.0, 5.2;
    el << 0.2, 0.9, 0.5, 1.2, 0.35;
    const gnss::Dop a = gnss::compute_dop(az, el);
    const gnss::Dop b = gnss::compute_dop((az.array() + 0.77).matrix(), el);
    CHECK(a.hdop == Approx(b.hdop));
    CHECK(a.vdop == Approx(b.vdop));
    CHECK(a.tdop == Approx(b.tdop));
}

TEST_CASE("DOP from az/el equals DOP from the rotated ECEF cofactor of chapter 1", "[dop][ch3]") {
    const geo::Geodetic site{geo::deg2rad(47.47), geo::deg2rad(19.06), 150.0};
    const Eigen::Vector3d r = geo::geodetic_to_ecef(site);
    Eigen::VectorXd az(7), el(7);
    az << 0.1, 0.9, 1.8, 2.5, 3.6, 4.4, 5.6;
    el << 0.25, 1.1, 0.6, 0.3, 0.8, 1.4, 0.45;
    // place_satellites uses the geocentric vertical; recompute the true az/el with the ellipsoidal ENU frame
    const Eigen::Matrix3Xd sats = gnss::place_satellites(r, az, el);
    const geo::LocalTangentPlane ltp(site);
    Eigen::VectorXd az_e(7), el_e(7);
    for (int i = 0; i < 7; ++i) {
        const Eigen::Vector3d d = ltp.ecef_to_enu(sats.col(i)).normalized();
        az_e(i) = std::atan2(d.x(), d.y());
        el_e(i) = std::asin(d.z());
    }
    const Eigen::MatrixX4d G = gnss::geometry_matrix(r, sats);
    const Eigen::Matrix4d Q = (G.transpose() * G).inverse();
    const gnss::Dop from_ecef = gnss::dop_from_ecef_cofactor(Q, site.lat, site.lon);
    const gnss::Dop from_azel = gnss::compute_dop(az_e, el_e);
    CHECK(from_ecef.hdop == Approx(from_azel.hdop).epsilon(1e-9));
    CHECK(from_ecef.vdop == Approx(from_azel.vdop).epsilon(1e-9));
    CHECK(from_ecef.gdop == Approx(from_azel.gdop).epsilon(1e-9));
    // trace (and so GDOP, PDOP) does not depend on the frame at all
    CHECK(std::sqrt(Q.trace()) == Approx(from_azel.gdop).epsilon(1e-9));
}

TEST_CASE("degenerate geometries are flagged", "[dop][ch3]") {
    Eigen::VectorXd az(3), el(3);
    az << 0.0, 2.0, 4.0;
    el << 0.5, 0.5, 0.5;
    CHECK_FALSE(gnss::compute_dop(az, el).valid);  // three satellites
    Eigen::VectorXd az4 = Eigen::VectorXd::LinSpaced(4, 0.0, 4.5), el4 = Eigen::VectorXd::Constant(4, 0.0);
    // all on the horizon: up and clock columns are identical up to sign -> singular
    CHECK_FALSE(gnss::compute_dop(az4, el4).valid);
}

TEST_CASE("elevation weights reduce to unit weights for equal elevations", "[dop][ch3]") {
    Eigen::VectorXd az(6), el = Eigen::VectorXd::Constant(6, 0.7);
    for (int i = 0; i < 6; ++i) az(i) = i;
    el(5) = kPi / 2;
    const Eigen::VectorXd w = gnss::elevation_weights(el);
    CHECK(w(5) == Approx(1.0));
    CHECK(w(0) == Approx(std::sin(0.7) * std::sin(0.7)));
    // scaling all weights by c scales Q by 1/c
    const gnss::Dop a = gnss::compute_dop(az, el, Eigen::VectorXd::Constant(6, 4.0));
    const gnss::Dop b = gnss::compute_dop(az, el);
    CHECK(a.hdop == Approx(b.hdop / 2.0));
}

TEST_CASE("error ellipse axes and probability content", "[ellipse][ch3]") {
    Eigen::Matrix2d cov;
    cov << 9.0, 0.0, 0.0, 4.0;
    const auto e1 = gnss::error_ellipse(cov);  // 1-sigma ellipse holds 39.35 %
    CHECK(e1.semi_major == Approx(3.0));
    CHECK(e1.semi_minor == Approx(2.0));
    CHECK(std::abs(std::sin(e1.angle)) < 1e-12);
    const auto e95 = gnss::error_ellipse(cov, 0.95);
    CHECK(e95.semi_major == Approx(3.0 * std::sqrt(5.991464547)).epsilon(1e-8));

    // rotated covariance: major axis along 30 degrees
    const double a = kPi / 6;
    Eigen::Matrix2d R;
    R << std::cos(a), -std::sin(a), std::sin(a), std::cos(a);
    const auto er = gnss::error_ellipse(R * cov * R.transpose());
    CHECK(std::abs(std::sin(er.angle - a)) < 1e-12);

    // Monte Carlo: the 95 % ellipse contains 95 % of samples
    std::mt19937 rng(5);
    std::normal_distribution<double> g(0.0, 1.0);
    const Eigen::Matrix2d L = (R * cov * R.transpose()).llt().matrixL();
    const auto e = gnss::error_ellipse(R * cov * R.transpose(), 0.95);
    int inside = 0, total = 20000;
    for (int i = 0; i < total; ++i) {
        const Eigen::Vector2d x = L * Eigen::Vector2d(g(rng), g(rng));
        const double u = std::cos(e.angle) * x.x() + std::sin(e.angle) * x.y();
        const double v = -std::sin(e.angle) * x.x() + std::cos(e.angle) * x.y();
        if (u * u / (e.semi_major * e.semi_major) + v * v / (e.semi_minor * e.semi_minor) <= 1.0) ++inside;
    }
    CHECK(inside / double(total) == Approx(0.95).margin(0.006));
}

TEST_CASE("precision statistics of an isotropic Gaussian cloud", "[precision][ch3]") {
    std::mt19937 rng(11);
    std::normal_distribution<double> g(0.0, 2.0);
    Eigen::Matrix2Xd en(2, 50000);
    for (int i = 0; i < en.cols(); ++i) en.col(i) << 10.0 + g(rng), -5.0 + g(rng);
    const auto s = gnss::precision_stats(en);
    CHECK(s.mean.x() == Approx(10.0).margin(0.05));
    CHECK(s.sigma_east == Approx(2.0).epsilon(0.02));
    CHECK(s.drms == Approx(2.0 * std::sqrt(2.0)).epsilon(0.02));
    // Rayleigh distribution: median sigma sqrt(2 ln 2), 95 % sigma sqrt(-2 ln 0.05)
    CHECK(s.cep50 == Approx(2.0 * std::sqrt(2 * std::log(2.0))).epsilon(0.02));
    CHECK(s.r95 == Approx(2.0 * std::sqrt(-2 * std::log(0.05))).epsilon(0.02));
}

TEST_CASE("motion classifier separates standing and walking", "[motion][ch3]") {
    // 60 s standing, 60 s walking east at 1.2 m/s, 60 s standing; 1 Hz fixes with 2 m noise
    std::mt19937 rng(2);
    std::normal_distribution<double> g(0.0, 2.0);
    const int n = 180;
    Eigen::VectorXd t(n), speed(n);
    Eigen::Matrix2Xd en(2, n);
    double x = 0.0;
    for (int i = 0; i < n; ++i) {
        t(i) = i;
        const bool walking = i >= 60 && i < 120;
        if (walking) x += 1.2;
        en.col(i) << x + g(rng), g(rng);
        speed(i) = std::numeric_limits<double>::quiet_NaN();  // position only
    }
    gnss::MotionParams p;
    p.sigma_displacement_m = 2.0 * std::sqrt(2.0);  // white noise: difference of two fixes
    p.window_s = 10.0;
    const auto r = gnss::classify_motion(t, en, speed, p);
    int errors = 0;
    for (int i = 15; i < n; ++i) {
        const bool truth = i >= 60 && i < 120;
        // the window delays detection by up to window_s at each transition
        if ((i >= 60 && i < 75) || (i >= 120 && i < 135)) continue;
        if (r.moving[static_cast<std::size_t>(i)] != truth) ++errors;
    }
    CHECK(errors == 0);
    CHECK(std::isnan(r.test_statistic(0)));
    CHECK(r.displacement_speed(100) == Approx(1.2).margin(0.6));

    // a receiver speed above the threshold alone marks motion
    Eigen::VectorXd fast = Eigen::VectorXd::Constant(n, 2.0);
    Eigen::Matrix2Xd still = Eigen::Matrix2Xd::Zero(2, n);
    const auto r2 = gnss::classify_motion(t, still, fast, p);
    CHECK(r2.moving[50]);
}

TEST_CASE("fix-code names", "[ch3]") {
    CHECK(gnss::gga_quality_name(4) == "RTK fixed");
    CHECK(gnss::gsa_fix_type_name(2) == "2D");
    CHECK(gnss::navsat_status_name(2) == "STATUS_GBAS_FIX");
}
