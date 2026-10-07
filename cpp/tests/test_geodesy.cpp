// Coordinate conversions against PROJ 9 (via pyproj 3.8) reference values, closed-form special cases and
// round trips (1_theory/02_coordinate_frames.md).
#include <Eigen/Dense>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <random>

#include "sensor_fusion/geodesy.hpp"

using Catch::Approx;
namespace geo = sensor_fusion::geo;
using geo::deg2rad;

namespace {
constexpr double kPi = 3.14159265358979323846;
geo::Geodetic deg(double lat, double lon, double h = 0.0) { return {deg2rad(lat), deg2rad(lon), h}; }
}  // namespace

TEST_CASE("WGS-84 derived constants", "[geo][ch2]") {
    CHECK(geo::kWgs84.b() == Approx(6356752.314245).margin(1e-6));
    CHECK(geo::kWgs84.e2() == Approx(0.00669437999014).epsilon(1e-12));
    CHECK(geo::prime_vertical_radius(0.0) == Approx(geo::kWgs84.a));
    CHECK(geo::meridian_radius(0.0) == Approx(geo::kWgs84.a * (1 - geo::kWgs84.e2())));
    // at the pole both radii equal a^2 / b
    const double polar = geo::kWgs84.a * geo::kWgs84.a / geo::kWgs84.b();
    CHECK(geo::prime_vertical_radius(deg2rad(90.0)) == Approx(polar));
    CHECK(geo::meridian_radius(deg2rad(90.0)) == Approx(polar));
}

TEST_CASE("geodetic -> ECEF: special points and PROJ reference values", "[geo][ch2]") {
    const auto eq = geo::geodetic_to_ecef(deg(0, 0, 0));
    CHECK(eq.isApprox(Eigen::Vector3d(geo::kWgs84.a, 0, 0)));
    const auto np = geo::geodetic_to_ecef(deg(90, 0, 0));
    CHECK(np.z() == Approx(geo::kWgs84.b()));
    CHECK(std::abs(np.x()) < 1e-9);

    struct Ref {
        double lat, lon, h, x, y, z;
    };
    const Ref refs[] = {
        {47.47439833333333, 19.056951666666667, 146.6, 4082357.096810, 1410208.773608, 4677682.316472},
        {37.7749, -122.4194, 15.0, -2706181.202951, -4261069.497794, 3885734.678436},
        {-33.8688, 151.2093, 40.0, -4646080.379115, 2553222.337804, -3534394.679635},
    };
    for (const auto& r : refs) {
        const auto p = geo::geodetic_to_ecef(deg(r.lat, r.lon, r.h));
        CHECK((p - Eigen::Vector3d(r.x, r.y, r.z)).norm() < 1e-5);
    }
}

TEST_CASE("ECEF -> geodetic round trip, surface to GNSS orbit, all latitudes", "[geo][ch2]") {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> lat(-90.0, 90.0), lon(-180.0, 180.0), h(-5000.0, 2.7e7);
    for (int i = 0; i < 2000; ++i) {
        const geo::Geodetic g = deg(lat(rng), lon(rng), h(rng));
        const geo::Geodetic back = geo::ecef_to_geodetic(geo::geodetic_to_ecef(g));
        CHECK(std::abs(back.lat - g.lat) < 1e-12);
        CHECK(std::abs(back.h - g.h) < 1e-5);
        if (std::abs(g.lat) < deg2rad(89.999)) CHECK(std::abs(std::remainder(back.lon - g.lon, 2 * kPi)) < 1e-12);
    }
    // exactly on the axis
    const auto pole = geo::ecef_to_geodetic(Eigen::Vector3d(0, 0, geo::kWgs84.b() + 10.0));
    CHECK(pole.lat == Approx(kPi / 2));
    CHECK(pole.h == Approx(10.0).margin(1e-6));
}

TEST_CASE("ENU rotation is orthonormal and points the right way", "[geo][ch2]") {
    const geo::Geodetic o = deg(47.47, 19.06, 150.0);
    const Eigen::Matrix3d R = geo::rotation_ecef_enu(o.lat, o.lon);
    CHECK((R.transpose() * R - Eigen::Matrix3d::Identity()).norm() < 1e-14);
    CHECK(R.determinant() == Approx(1.0));
    // "up" is the ellipsoid normal: moving up in h moves along the third column
    const geo::LocalTangentPlane ltp(o);
    CHECK(ltp.geodetic_to_enu(o).norm() < 1e-9);
    const auto up = ltp.geodetic_to_enu(deg(47.47, 19.06, 250.0));
    CHECK(up.isApprox(Eigen::Vector3d(0, 0, 100.0), 1e-9));
    // a small step north increases N only (to first order)
    const auto north = ltp.geodetic_to_enu(deg(47.4701, 19.06, 150.0));
    CHECK(north.y() > 11.0);
    CHECK(std::abs(north.x()) < 1e-6);
    // round trip
    const Eigen::Vector3d p(1234.5, -987.6, 55.5);
    CHECK((ltp.geodetic_to_enu(ltp.enu_to_geodetic(p)) - p).norm() < 1e-6);
}

TEST_CASE("flat-Earth approximation: exact to first order, error grows quadratically", "[geo][ch2]") {
    const geo::Geodetic o = deg(47.47, 19.06, 150.0);
    const geo::LocalTangentPlane ltp(o);
    auto error_at = [&](double d) {  // point d metres north-east of the origin
        const geo::Geodetic g = ltp.enu_to_geodetic(Eigen::Vector3d(d, d, 0.0));
        return (geo::geodetic_to_enu_flat(g, o) - ltp.geodetic_to_enu(g)).norm();
    };
    // second-order terms of eq. (2.7): up (e^2 + n^2) / 2R and east e n tan(phi0) / R
    const double R = 6.371e6, d = 100.0;
    const double predicted = std::hypot((d * d + d * d) / (2 * R), d * d * std::tan(o.lat) / R);
    CHECK(error_at(d) == Approx(predicted).epsilon(0.1));
    const double e1 = error_at(1000.0), e10 = error_at(10000.0);
    CHECK(e10 / e1 == Approx(100.0).epsilon(0.1));  // quadratic in distance
}

TEST_CASE("UTM zones, including Norway and Svalbard", "[geo][ch2]") {
    CHECK(geo::utm_zone(deg2rad(47.47), deg2rad(19.06)) == 34);
    CHECK(geo::utm_zone(deg2rad(0.0), deg2rad(-180.0)) == 1);
    CHECK(geo::utm_zone(deg2rad(0.0), deg2rad(179.9)) == 60);
    CHECK(geo::utm_zone(deg2rad(60.0), deg2rad(5.0)) == 32);   // Bergen: 32V, not 31V
    CHECK(geo::utm_zone(deg2rad(78.0), deg2rad(15.0)) == 33);  // Longyearbyen
    CHECK(geo::utm_zone(deg2rad(78.0), deg2rad(10.0)) == 33);
    CHECK(geo::utm_zone(deg2rad(78.0), deg2rad(8.0)) == 31);
    CHECK(geo::utm_central_meridian(34) == Approx(deg2rad(21.0)));
}

TEST_CASE("geodetic -> UTM against PROJ, including convergence and scale", "[geo][ch2]") {
    struct Ref {
        double lat, lon;
        int zone;
        double e, n;
    };
    const Ref refs[] = {
        {47.47439833333333, 19.056951666666667, 34, 353593.626571, 5259714.355562},
        {37.7749, -122.4194, 10, 551130.768481, 4180998.881499},
        {-33.8688, 151.2093, 56, 334368.633648, 6250948.345385},
        {0.0, 3.0, 31, 500000.0, 0.0},
        {51.5074, -0.1278, 30, 699316.234312, 5710163.758081},
        {-33.4489, -70.6693, 19, 344846.720310, 6297700.155610},
    };
    for (const auto& r : refs) {
        const geo::Utm u = geo::geodetic_to_utm(deg2rad(r.lat), deg2rad(r.lon));
        INFO("lat " << r.lat << " lon " << r.lon);
        CHECK(u.zone == r.zone);
        CHECK(u.north == (r.lat >= 0));
        CHECK(std::abs(u.easting - r.e) < 1e-5);
        CHECK(std::abs(u.northing - r.n) < 1e-5);
        const geo::Geodetic back = geo::utm_to_geodetic(u);
        CHECK(std::abs(back.lat - deg2rad(r.lat)) < 1e-13);
        CHECK(std::abs(back.lon - deg2rad(r.lon)) < 1e-13);
    }
    const geo::Utm b = geo::geodetic_to_utm(deg2rad(47.47439833333333), deg2rad(19.056951666666667));
    CHECK(geo::rad2deg(b.convergence) == Approx(-1.4322319399).margin(1e-9));
    CHECK(b.scale == Approx(0.999863418866).margin(1e-11));
    const geo::Utm s = geo::geodetic_to_utm(deg2rad(-33.8688), deg2rad(151.2093));
    CHECK(geo::rad2deg(s.convergence) == Approx(0.9981718559).margin(1e-9));
    CHECK(s.scale == Approx(0.999938200544).margin(1e-11));
    // on the central meridian the scale is k0 and the convergence zero
    const geo::Utm c = geo::geodetic_to_utm(deg2rad(30.0), deg2rad(21.0));
    CHECK(c.scale == Approx(0.9996).margin(1e-12));
    CHECK(std::abs(c.convergence) < 1e-15);
    CHECK(c.easting == Approx(500000.0).margin(1e-9));
}

TEST_CASE("UTM round trip across a zone, both hemispheres", "[geo][ch2]") {
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> lat(-80.0, 84.0), dlon(-3.0, 3.0);
    for (int i = 0; i < 2000; ++i) {
        const double la = deg2rad(lat(rng));
        const double lo = deg2rad(21.0 + dlon(rng));
        const geo::Utm u = geo::geodetic_to_utm(la, lo, 34);
        const geo::Geodetic g = geo::utm_to_geodetic(u);
        CHECK(std::abs(g.lat - la) < 1e-12);
        CHECK(std::abs(g.lon - lo) < 1e-12);
    }
}
