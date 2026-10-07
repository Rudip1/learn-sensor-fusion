#include "sensor_fusion/geodesy.hpp"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <stdexcept>

namespace sensor_fusion::geo {

namespace {

constexpr double kPi = 3.14159265358979323846;

/// Krueger series coefficients for a given third flattening n (Karney 2011, eqs. 35 and 36).
struct KruegerSeries {
    double A;                    ///< rectifying radius over a, times a: A = a/(1+n) (1 + n^2/4 + ...)
    std::array<double, 6> alpha;  ///< forward
    std::array<double, 6> beta;   ///< inverse
};

KruegerSeries krueger(const Ellipsoid& ell) {
    const double n = ell.f / (2.0 - ell.f);
    const double n2 = n * n, n3 = n2 * n, n4 = n3 * n, n5 = n4 * n, n6 = n5 * n;
    KruegerSeries s;
    s.A = ell.a / (1.0 + n) * (1.0 + n2 / 4.0 + n4 / 64.0 + n6 / 256.0);  // eq. (2.9)
    s.alpha = {n / 2 - 2 * n2 / 3 + 5 * n3 / 16 + 41 * n4 / 180 - 127 * n5 / 288 + 7891 * n6 / 37800,
               13 * n2 / 48 - 3 * n3 / 5 + 557 * n4 / 1440 + 281 * n5 / 630 - 1983433 * n6 / 1935360,
               61 * n3 / 240 - 103 * n4 / 140 + 15061 * n5 / 26880 + 167603 * n6 / 181440,
               49561 * n4 / 161280 - 179 * n5 / 168 + 6601661 * n6 / 7257600,
               34729 * n5 / 80640 - 3418889 * n6 / 1995840,
               212378941 * n6 / 319334400};
    s.beta = {n / 2 - 2 * n2 / 3 + 37 * n3 / 96 - n4 / 360 - 81 * n5 / 512 + 96199 * n6 / 604800,
              n2 / 48 + n3 / 15 - 437 * n4 / 1440 + 46 * n5 / 105 - 1118711 * n6 / 3870720,
              17 * n3 / 480 - 37 * n4 / 840 - 209 * n5 / 4480 + 5569 * n6 / 90720,
              4397 * n4 / 161280 - 11 * n5 / 504 - 830251 * n6 / 7257600,
              4583 * n5 / 161280 - 108847 * n6 / 3991680,
              20648693 * n6 / 638668800};
    return s;
}

/// tan of the conformal latitude from tan of the geodetic latitude. Eq. (2.8).
double conformal_tan(double tau, double e) {
    const double sigma = std::sinh(e * std::atanh(e * tau / std::hypot(1.0, tau)));
    return tau * std::hypot(1.0, sigma) - sigma * std::hypot(1.0, tau);
}

/// Inverse of conformal_tan by Newton's method (Karney 2011, eqs. 19-21).
double geodetic_tan(double tau_prime, const Ellipsoid& ell) {
    const double e2 = ell.e2(), e = std::sqrt(e2);
    double tau = tau_prime / (1.0 - e2);  // good start: exact for the polar limit
    for (int i = 0; i < 8; ++i) {
        const double tp = conformal_tan(tau, e);
        const double dtau = (tau_prime - tp) * (1.0 + (1.0 - e2) * tau * tau) /
                            ((1.0 - e2) * std::hypot(1.0, tp) * std::hypot(1.0, tau));
        tau += dtau;
        if (std::abs(dtau) < 1e-14 * std::max(1.0, std::abs(tau))) break;
    }
    return tau;
}

}  // namespace

double prime_vertical_radius(double lat, const Ellipsoid& ell) {
    const double s = std::sin(lat);
    return ell.a / std::sqrt(1.0 - ell.e2() * s * s);  // eq. (2.2)
}

double meridian_radius(double lat, const Ellipsoid& ell) {
    const double s = std::sin(lat);
    const double w2 = 1.0 - ell.e2() * s * s;
    return ell.a * (1.0 - ell.e2()) / (w2 * std::sqrt(w2));  // eq. (2.6)
}

Eigen::Vector3d geodetic_to_ecef(const Geodetic& g, const Ellipsoid& ell) {
    const double N = prime_vertical_radius(g.lat, ell);
    const double cl = std::cos(g.lat), sl = std::sin(g.lat);
    return {(N + g.h) * cl * std::cos(g.lon),  // eq. (2.3)
            (N + g.h) * cl * std::sin(g.lon),
            (N * (1.0 - ell.e2()) + g.h) * sl};
}

Geodetic ecef_to_geodetic(const Eigen::Vector3d& p, const Ellipsoid& ell) {
    const double e2 = ell.e2();
    const double rho = std::hypot(p.x(), p.y());  // distance from the polar axis
    Geodetic g;
    g.lon = std::atan2(p.y(), p.x());
    // Algorithm 2.1: start from the geocentric-like guess, iterate phi = atan2(z + e^2 N sin(phi), rho).
    double lat = std::atan2(p.z(), rho * (1.0 - e2));
    for (int i = 0; i < 30; ++i) {
        const double N = prime_vertical_radius(lat, ell);
        const double next = std::atan2(p.z() + e2 * N * std::sin(lat), rho);
        const bool done = std::abs(next - lat) < 1e-15;
        lat = next;
        if (done) break;
    }
    g.lat = lat;
    // h = rho cos(phi) + z sin(phi) - a sqrt(1 - e^2 sin^2 phi): exact and well conditioned at all latitudes.
    const double s = std::sin(lat);
    g.h = rho * std::cos(lat) + p.z() * s - ell.a * std::sqrt(1.0 - e2 * s * s);
    return g;
}

Eigen::Matrix3d rotation_ecef_enu(double lat, double lon) {
    const double sl = std::sin(lat), cl = std::cos(lat), so = std::sin(lon), co = std::cos(lon);
    Eigen::Matrix3d R;
    // columns: east, north, up (eq. 2.4)
    R << -so, -sl * co, cl * co,
          co, -sl * so, cl * so,
         0.0,       cl,      sl;
    return R;
}

LocalTangentPlane::LocalTangentPlane(const Geodetic& origin, const Ellipsoid& ell)
    : ell_(ell), origin_(origin), origin_ecef_(geodetic_to_ecef(origin, ell)),
      R_(rotation_ecef_enu(origin.lat, origin.lon)) {}

Eigen::Vector3d LocalTangentPlane::ecef_to_enu(const Eigen::Vector3d& p_ecef) const {
    return R_.transpose() * (p_ecef - origin_ecef_);  // eq. (2.5)
}

Eigen::Vector3d LocalTangentPlane::enu_to_ecef(const Eigen::Vector3d& p_enu) const {
    return origin_ecef_ + R_ * p_enu;
}

Eigen::Vector3d LocalTangentPlane::geodetic_to_enu(const Geodetic& g) const {
    return ecef_to_enu(geodetic_to_ecef(g, ell_));
}

Geodetic LocalTangentPlane::enu_to_geodetic(const Eigen::Vector3d& p_enu) const {
    return ecef_to_geodetic(enu_to_ecef(p_enu), ell_);
}

Eigen::Vector3d geodetic_to_enu_flat(const Geodetic& g, const Geodetic& origin, const Ellipsoid& ell) {
    double dlon = g.lon - origin.lon;
    if (dlon > kPi) dlon -= 2 * kPi;
    if (dlon < -kPi) dlon += 2 * kPi;
    const double N = prime_vertical_radius(origin.lat, ell);
    const double M = meridian_radius(origin.lat, ell);
    return {(N + origin.h) * std::cos(origin.lat) * dlon,  // eq. (2.7)
            (M + origin.h) * (g.lat - origin.lat),
            g.h - origin.h};
}

int utm_zone(double lat, double lon) {
    const double lat_d = rad2deg(lat);
    double lon_d = rad2deg(lon);
    lon_d = lon_d - 360.0 * std::floor((lon_d + 180.0) / 360.0);  // wrap to [-180, 180)
    int zone = static_cast<int>(std::floor((lon_d + 180.0) / 6.0)) + 1;
    if (zone > 60) zone = 60;
    if (lat_d >= 56.0 && lat_d < 64.0 && lon_d >= 3.0 && lon_d < 12.0) zone = 32;  // south-west Norway
    if (lat_d >= 72.0 && lat_d < 84.0 && lon_d >= 0.0 && lon_d < 42.0) {            // Svalbard
        if (lon_d < 9.0)
            zone = 31;
        else if (lon_d < 21.0)
            zone = 33;
        else if (lon_d < 33.0)
            zone = 35;
        else
            zone = 37;
    }
    return zone;
}

double utm_central_meridian(int zone) {
    if (zone < 1 || zone > 60) throw std::invalid_argument("UTM zone must be in 1..60");
    return deg2rad(6.0 * zone - 183.0);
}

Eigen::Vector4d TransverseMercator::forward(double lat, double lon) const {
    const KruegerSeries s = krueger(ell);
    const double e = std::sqrt(ell.e2());
    double dlon = lon - central_meridian;
    dlon = std::remainder(dlon, 2 * kPi);

    // Conformal latitude, then the spherical transverse Mercator (Gauss-Schreiber). Eq. (2.8).
    const double tau = std::tan(lat);
    const double tau_p = conformal_tan(tau, e);
    const double cl = std::cos(dlon), sl = std::sin(dlon);
    const double xi_p = std::atan2(tau_p, cl);
    const double eta_p = std::asinh(sl / std::hypot(tau_p, cl));

    // Krueger series, eq. (2.10), and its derivative for convergence and scale.
    double xi = xi_p, eta = eta_p, p = 1.0, q = 0.0;
    for (int j = 1; j <= 6; ++j) {
        const double a = s.alpha[j - 1];
        const double c2 = std::cos(2 * j * xi_p), s2 = std::sin(2 * j * xi_p);
        const double ch = std::cosh(2 * j * eta_p), sh = std::sinh(2 * j * eta_p);
        xi += a * s2 * ch;
        eta += a * c2 * sh;
        p += 2 * j * a * c2 * ch;
        q += 2 * j * a * s2 * sh;
    }
    const double easting = false_easting + k0 * s.A * eta;  // eq. (2.11)
    const double northing = false_northing + k0 * s.A * xi;

    // Grid convergence and point scale (Karney 2011, eqs. 13-14 and 22-25).
    const double gamma = std::atan(tau_p / std::hypot(1.0, tau_p) * std::tan(dlon)) + std::atan2(q, p);
    const double sphi = std::sin(lat);
    const double k = k0 * std::sqrt(1.0 - ell.e2() * sphi * sphi) * std::hypot(1.0, tau) /
                     std::hypot(tau_p, cl) * (s.A / ell.a) * std::hypot(p, q);
    return {easting, northing, gamma, k};
}

Eigen::Vector2d TransverseMercator::inverse(double easting, double northing) const {
    const KruegerSeries s = krueger(ell);
    const double xi = (northing - false_northing) / (k0 * s.A);
    const double eta = (easting - false_easting) / (k0 * s.A);
    double xi_p = xi, eta_p = eta;
    for (int j = 1; j <= 6; ++j) {  // eq. (2.12)
        const double b = s.beta[j - 1];
        xi_p -= b * std::sin(2 * j * xi) * std::cosh(2 * j * eta);
        eta_p -= b * std::cos(2 * j * xi) * std::sinh(2 * j * eta);
    }
    const double sinh_eta = std::sinh(eta_p), cos_xi = std::cos(xi_p);
    const double tau_p = std::sin(xi_p) / std::hypot(sinh_eta, cos_xi);
    const double dlon = std::atan2(sinh_eta, cos_xi);
    const double lat = std::atan(geodetic_tan(tau_p, ell));
    return {lat, central_meridian + dlon};
}

Utm geodetic_to_utm(double lat, double lon, int zone, const Ellipsoid& ell) {
    Utm u;
    u.zone = zone > 0 ? zone : utm_zone(lat, lon);
    u.north = lat >= 0.0;
    TransverseMercator tm;
    tm.central_meridian = utm_central_meridian(u.zone);
    tm.false_northing = u.north ? 0.0 : 10'000'000.0;
    tm.ell = ell;
    const Eigen::Vector4d r = tm.forward(lat, lon);
    u.easting = r(0);
    u.northing = r(1);
    u.convergence = r(2);
    u.scale = r(3);
    return u;
}

Geodetic utm_to_geodetic(const Utm& u, const Ellipsoid& ell) {
    TransverseMercator tm;
    tm.central_meridian = utm_central_meridian(u.zone);
    tm.false_northing = u.north ? 0.0 : 10'000'000.0;
    tm.ell = ell;
    const Eigen::Vector2d ll = tm.inverse(u.easting, u.northing);
    Geodetic g;
    g.lat = ll(0);
    g.lon = std::remainder(ll(1), 2 * kPi);
    return g;
}

}  // namespace sensor_fusion::geo
