#pragma once
/// @file geodesy.hpp
/// WGS-84 ellipsoid, geodetic <-> ECEF <-> local ENU conversions and the UTM projection.
/// Theory: 1_theory/02_coordinate_frames.md.

#include <Eigen/Core>
#include <string>

namespace sensor_fusion::geo {

/// Reference ellipsoid, defined by semi-major axis a and flattening f. Eq. (2.1).
struct Ellipsoid {
    double a;  ///< semi-major axis [m]
    double f;  ///< flattening

    constexpr double b() const { return a * (1.0 - f); }       ///< semi-minor axis
    constexpr double e2() const { return f * (2.0 - f); }      ///< first eccentricity squared
    constexpr double ep2() const { return e2() / (1.0 - e2()); }  ///< second eccentricity squared
};

/// WGS-84, the ellipsoid of GPS (NIMA TR8350.2).
inline constexpr Ellipsoid kWgs84{6378137.0, 1.0 / 298.257223563};

/// Geodetic coordinates: latitude and longitude in radians, ellipsoidal height in metres.
struct Geodetic {
    double lat = 0.0;
    double lon = 0.0;
    double h = 0.0;
};

/// Degrees <-> radians.
constexpr double deg2rad(double deg) { return deg * 0.017453292519943295; }
constexpr double rad2deg(double rad) { return rad * 57.29577951308232; }

/// Prime-vertical radius of curvature N(phi). Eq. (2.2).
double prime_vertical_radius(double lat, const Ellipsoid& ell = kWgs84);
/// Meridian radius of curvature M(phi). Eq. (2.6).
double meridian_radius(double lat, const Ellipsoid& ell = kWgs84);

/// (phi, lambda, h) -> ECEF. Eq. (2.3).
Eigen::Vector3d geodetic_to_ecef(const Geodetic& g, const Ellipsoid& ell = kWgs84);

/// ECEF -> (phi, lambda, h) by fixed-point iteration on the latitude (Algorithm 2.1). Converges to
/// below 1e-12 rad in a few iterations for points from the Earth's centre region to GNSS orbits.
Geodetic ecef_to_geodetic(const Eigen::Vector3d& p, const Ellipsoid& ell = kWgs84);

/// Rotation R_EN whose columns are the east, north and up unit vectors at (lat, lon), expressed in ECEF.
/// ECEF -> ENU is R_EN^T. Eq. (2.4).
Eigen::Matrix3d rotation_ecef_enu(double lat, double lon);

/// A local east-north-up frame {N} with its origin at a geodetic point. Eq. (2.5).
class LocalTangentPlane {
public:
    explicit LocalTangentPlane(const Geodetic& origin, const Ellipsoid& ell = kWgs84);

    Eigen::Vector3d ecef_to_enu(const Eigen::Vector3d& p_ecef) const;
    Eigen::Vector3d enu_to_ecef(const Eigen::Vector3d& p_enu) const;
    Eigen::Vector3d geodetic_to_enu(const Geodetic& g) const;
    Geodetic enu_to_geodetic(const Eigen::Vector3d& p_enu) const;

    const Geodetic& origin() const { return origin_; }
    const Eigen::Vector3d& origin_ecef() const { return origin_ecef_; }
    const Eigen::Matrix3d& rotation() const { return R_; }  ///< R_EN

private:
    Ellipsoid ell_;
    Geodetic origin_;
    Eigen::Vector3d origin_ecef_;
    Eigen::Matrix3d R_;
};

/// First-order ("flat Earth") local coordinates: east = N(phi0) cos(phi0) dlambda, north = M(phi0) dphi,
/// up = dh. Eq. (2.7). Valid for small areas; the error grows with the square of the distance.
Eigen::Vector3d geodetic_to_enu_flat(const Geodetic& g, const Geodetic& origin, const Ellipsoid& ell = kWgs84);

/// A position in the Universal Transverse Mercator grid.
struct Utm {
    int zone = 0;          ///< 1..60
    bool north = true;     ///< hemisphere (false: false northing of 10 000 km)
    double easting = 0.0;  ///< [m], false easting 500 km at the central meridian
    double northing = 0.0; ///< [m]
    double convergence = 0.0;  ///< grid convergence gamma [rad]: angle from true north to grid north
    double scale = 0.0;        ///< point scale factor k
};

/// Standard UTM zone of a point, including the Norway (32V) and Svalbard (31X-37X) exceptions.
int utm_zone(double lat, double lon);

/// Central meridian of a UTM zone [rad].
double utm_central_meridian(int zone);

/// Geodetic -> UTM with the Krueger series to sixth order in n (Karney 2011). Eqs. (2.8)-(2.11).
/// @param zone force a zone (0 = standard zone of the point).
Utm geodetic_to_utm(double lat, double lon, int zone = 0, const Ellipsoid& ell = kWgs84);

/// UTM -> (lat, lon) [rad]. Inverse series (2.12) and Newton iteration for the latitude.
Geodetic utm_to_geodetic(const Utm& u, const Ellipsoid& ell = kWgs84);

/// Transverse Mercator with arbitrary central meridian, scale and false origin (UTM is a special case).
struct TransverseMercator {
    double central_meridian = 0.0;  ///< [rad]
    double k0 = 0.9996;
    double false_easting = 500000.0;
    double false_northing = 0.0;
    Ellipsoid ell = kWgs84;

    /// (lat, lon) [rad] -> (E, N, gamma, k).
    Eigen::Vector4d forward(double lat, double lon) const;
    /// (E, N) -> (lat, lon) [rad].
    Eigen::Vector2d inverse(double easting, double northing) const;
};

}  // namespace sensor_fusion::geo
