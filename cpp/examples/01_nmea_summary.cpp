// Decode an NMEA 0183 log and print one line per epoch plus a summary.
//
//   01_nmea_summary data/gnss/nmea_log2.nmea [max_lines]
#include <cstdio>
#include <cstdlib>
#include <exception>

#include "sensor_fusion/nmea.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s LOG.nmea [max_lines]\n", argv[0]);
        return 1;
    }
    const int max_lines = argc > 2 ? std::atoi(argv[2]) : 10;
    try {
        const auto log = sensor_fusion::nmea::parse_log_file(argv[1]);
        std::printf("%-12s %4s %4s %4s %13s %13s %8s %6s %6s\n", "utc[s]", "qual", "used", "view", "lat[deg]",
                    "lon[deg]", "H[m]", "HDOP", "v[m/s]");
        int printed = 0;
        for (const auto& e : log.epochs) {
            if (printed++ >= max_lines) break;
            std::printf("%-12.3f %4d %4d %4zu %13.7f %13.7f %8.1f %6.2f %6.2f\n", e.time_of_day_s.value_or(-1.0),
                        e.fix_quality, e.satellites_used, e.satellites.size(), e.latitude_deg.value_or(0.0),
                        e.longitude_deg.value_or(0.0), e.altitude_msl_m.value_or(0.0), e.hdop.value_or(0.0),
                        e.speed_mps.value_or(0.0));
        }
        std::printf("\nlines %d, sentences %d, checksum errors %d, other errors %d, epochs %zu\n", log.stats.lines,
                    log.stats.sentences, log.stats.checksum_errors, log.stats.other_errors, log.epochs.size());
        for (const auto& [type, count] : log.stats.counts_by_type) std::printf("  %s: %d\n", type.c_str(), count);
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "error: %s\n", ex.what());
        return 1;
    }
    return 0;
}
