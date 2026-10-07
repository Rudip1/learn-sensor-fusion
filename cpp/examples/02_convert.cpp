// Convert a geodetic position to ECEF, to ENU about an origin, and to UTM.
//
//   02_convert LAT_DEG LON_DEG H_M [LAT0_DEG LON0_DEG H0_M]
#include <cstdio>
#include <cstdlib>

#include "sensor_fusion/geodesy.hpp"

int main(int argc, char** argv) {
    namespace geo = sensor_fusion::geo;
    if (argc != 4 && argc != 7) {
        std::fprintf(stderr, "usage: %s LAT_DEG LON_DEG H_M [LAT0_DEG LON0_DEG H0_M]\n", argv[0]);
        return 1;
    }
    const geo::Geodetic p{geo::deg2rad(std::atof(argv[1])), geo::deg2rad(std::atof(argv[2])), std::atof(argv[3])};
    const Eigen::Vector3d ecef = geo::geodetic_to_ecef(p);
    std::printf("ECEF  x %.4f  y %.4f  z %.4f m\n", ecef.x(), ecef.y(), ecef.z());

    const geo::Geodetic back = geo::ecef_to_geodetic(ecef);
    std::printf("back  lat %.10f  lon %.10f  h %.4f\n", geo::rad2deg(back.lat), geo::rad2deg(back.lon), back.h);

    const geo::Utm u = geo::geodetic_to_utm(p.lat, p.lon);
    std::printf("UTM   zone %d%c  E %.4f  N %.4f m   k %.7f  gamma %.5f deg\n", u.zone, u.north ? 'N' : 'S',
                u.easting, u.northing, u.scale, geo::rad2deg(u.convergence));

    if (argc == 7) {
        const geo::LocalTangentPlane ltp(
            {geo::deg2rad(std::atof(argv[4])), geo::deg2rad(std::atof(argv[5])), std::atof(argv[6])});
        const Eigen::Vector3d enu = ltp.geodetic_to_enu(p);
        const Eigen::Vector3d flat = geo::geodetic_to_enu_flat(p, ltp.origin());
        std::printf("ENU   e %.4f  n %.4f  u %.4f m\n", enu.x(), enu.y(), enu.z());
        std::printf("flat  e %.4f  n %.4f  u %.4f m   (difference %.4f m)\n", flat.x(), flat.y(), flat.z(),
                    (flat - enu).norm());
    }
    return 0;
}
