// Recompute DOP from the GSV azimuths/elevations of the satellites listed in GSA and compare it with the
// DOP the receiver reported, epoch by epoch.
//
//   03_dop_check data/gnss/nmea_log2.nmea
#include <cmath>
#include <cstdio>
#include <exception>
#include <map>

#include "sensor_fusion/gnss_quality.hpp"
#include "sensor_fusion/nmea.hpp"

int main(int argc, char** argv) {
    namespace nmea = sensor_fusion::nmea;
    namespace gnss = sensor_fusion::gnss;
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s LOG.nmea\n", argv[0]);
        return 1;
    }
    try {
        const auto log = nmea::parse_log_file(argv[1]);
        std::map<int, nmea::SatelliteInView> sky;  // latest GSV entry per PRN
        int compared = 0, agree = 0;
        std::printf("%10s %5s %13s %13s %13s\n", "utc[s]", "used", "PDOP rep/calc", "HDOP rep/calc", "VDOP rep/calc");
        for (const auto& e : log.epochs) {
            if (!e.satellites.empty()) {
                sky.clear();
                for (const auto& s : e.satellites)
                    if (s.azimuth_deg && s.elevation_deg) sky[s.prn] = s;
            }
            Eigen::VectorXd az(static_cast<Eigen::Index>(e.used_prns.size())), el(az.size());
            Eigen::Index n = 0;
            for (int prn : e.used_prns) {
                const auto it = sky.find(prn);
                if (it == sky.end()) break;
                az(n) = *it->second.azimuth_deg * 0.017453292519943295;
                el(n) = *it->second.elevation_deg * 0.017453292519943295;
                ++n;
            }
            if (n < 4 || n != az.size() || !e.pdop) continue;
            const gnss::Dop d = gnss::compute_dop(az, el);
            ++compared;
            agree += std::abs(d.pdop - *e.pdop) < 0.1;
            std::printf("%10.1f %5ld %6.2f/%6.2f %6.2f/%6.2f %6.2f/%6.2f\n", e.time_of_day_s.value_or(-1.0),
                        static_cast<long>(n), *e.pdop, d.pdop, e.hdop.value_or(NAN), d.hdop, e.vdop.value_or(NAN),
                        d.vdop);
        }
        std::printf("\nreported PDOP within 0.1 of the recomputed value in %d of %d epochs\n", agree, compared);
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "error: %s\n", ex.what());
        return 1;
    }
    return 0;
}
