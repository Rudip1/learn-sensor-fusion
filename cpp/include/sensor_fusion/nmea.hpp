#pragma once
/// @file nmea.hpp
/// NMEA 0183 sentence splitting, checksum and decoding of GGA, RMC, GSA, GSV and VTG.
/// Theory: 1_theory/01_gnss_fundamentals.md, section "NMEA 0183".

#include <cstdint>
#include <istream>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace sensor_fusion::nmea {

/// Outcome of splitting one line into an NMEA sentence.
enum class Status {
    ok,
    empty,               ///< blank line
    no_start_delimiter,  ///< no '$' in the line
    missing_checksum,    ///< no '*hh' suffix although one was required
    bad_checksum_format, ///< '*' not followed by exactly two hex digits
    checksum_mismatch,   ///< transmitted and computed checksum differ
    malformed            ///< address field too short to hold talker + type
};

/// Human-readable name of a status value.
const char* to_string(Status status);

/// One sentence split at the commas: "$GPGGA,a,b*hh" -> talker "GP", type "GGA", fields {"a","b"}.
struct Sentence {
    std::string talker;               ///< two-letter talker ID: GP, GL, GA, GB, GN, ...
    std::string type;                 ///< three-letter sentence formatter: GGA, RMC, ...
    std::vector<std::string> fields;  ///< data fields after the address field (may be empty strings)
};

struct SplitResult {
    Status status = Status::empty;
    Sentence sentence;
};

/// XOR of all bytes in @p body, the characters strictly between '$' and '*'. Eq. (1.6).
std::uint8_t checksum(std::string_view body);

/// Split a line into a sentence and verify its checksum. Leading and trailing whitespace are ignored; the
/// sentence starts at the first '$'.
SplitResult split(std::string_view line, bool require_checksum = true);

/// "ddmm.mmmm" / "dddmm.mmmm" plus hemisphere 'N','S','E','W' -> signed decimal degrees. Eq. (1.7).
/// Returns nullopt when either field is empty or malformed.
std::optional<double> parse_angle(std::string_view value, std::string_view hemisphere);

/// "hhmmss.sss" -> seconds since UTC midnight.
std::optional<double> parse_time_of_day(std::string_view value);

struct Date {
    int year = 0;  ///< four-digit year (two-digit NMEA years are mapped to 2000-2099)
    int month = 0;
    int day = 0;
};

/// "ddmmyy" -> date.
std::optional<Date> parse_date(std::string_view value);

/// Knots to metres per second, eq. (1.9).
constexpr double knots_to_mps(double knots) { return knots * 1852.0 / 3600.0; }

/// GGA: time, position and fix data.
struct Gga {
    std::optional<double> time_of_day_s;
    std::optional<double> latitude_deg;
    std::optional<double> longitude_deg;
    int fix_quality = 0;  ///< 0 invalid, 1 GPS (SPS), 2 DGPS, 4 RTK fixed, 5 RTK float, 6 dead reckoning
    int satellites_used = 0;
    std::optional<double> hdop;
    std::optional<double> altitude_msl_m;      ///< orthometric height H above mean sea level (geoid)
    std::optional<double> geoid_separation_m;  ///< N, height of the geoid above the ellipsoid
    std::optional<double> dgps_age_s;
    std::optional<int> dgps_station;
};

/// RMC: recommended minimum data.
struct Rmc {
    std::optional<double> time_of_day_s;
    bool valid = false;  ///< status 'A' (valid) vs 'V' (void)
    std::optional<double> latitude_deg;
    std::optional<double> longitude_deg;
    std::optional<double> speed_mps;   ///< speed over ground, converted from knots
    std::optional<double> course_deg;  ///< course over ground, degrees clockwise from true north
    std::optional<Date> date;
    char mode = 0;  ///< NMEA 2.3+: A autonomous, D differential, E estimated, N not valid
};

/// GSA: DOP and the satellites used in the solution.
struct Gsa {
    char selection_mode = 0;  ///< 'A' automatic 2D/3D, 'M' manual
    int fix_type = 1;         ///< 1 no fix, 2 2D, 3 3D
    std::vector<int> prns;    ///< PRNs of the satellites used (empty slots dropped)
    std::optional<double> pdop;
    std::optional<double> hdop;
    std::optional<double> vdop;
};

struct SatelliteInView {
    int prn = 0;
    std::optional<double> elevation_deg;
    std::optional<double> azimuth_deg;
    std::optional<double> snr_dbhz;  ///< carrier-to-noise density C/N0, empty when not tracked
};

/// GSV: satellites in view, up to four per sentence.
struct Gsv {
    int total_messages = 0;
    int message_number = 0;
    int satellites_in_view = 0;
    std::vector<SatelliteInView> satellites;
};

/// VTG: course and speed over ground.
struct Vtg {
    std::optional<double> course_true_deg;
    std::optional<double> course_magnetic_deg;
    std::optional<double> speed_mps;  ///< from the km/h field if present, else from knots
    char mode = 0;
};

std::optional<Gga> parse_gga(const Sentence& s);
std::optional<Rmc> parse_rmc(const Sentence& s);
std::optional<Gsa> parse_gsa(const Sentence& s);
std::optional<Gsv> parse_gsv(const Sentence& s);
std::optional<Vtg> parse_vtg(const Sentence& s);

/// Everything the receiver reported for one measurement epoch (one time tag).
struct Epoch {
    std::optional<double> time_of_day_s;
    std::optional<Date> date;
    // GGA
    int fix_quality = 0;
    int satellites_used = 0;
    std::optional<double> latitude_deg;
    std::optional<double> longitude_deg;
    std::optional<double> altitude_msl_m;
    std::optional<double> geoid_separation_m;
    std::optional<double> hdop;
    // GSA
    int fix_type = 1;
    std::vector<int> used_prns;
    std::optional<double> pdop;
    std::optional<double> vdop;
    // RMC / VTG
    bool rmc_valid = false;
    std::optional<double> speed_mps;
    std::optional<double> course_deg;
    // GSV
    std::vector<SatelliteInView> satellites;
};

struct LogStats {
    int lines = 0;
    int sentences = 0;  ///< lines that split with a valid checksum
    int checksum_errors = 0;
    int other_errors = 0;  ///< empty address, missing '$', bad fields
    std::map<std::string, int> counts_by_type;
};

struct Log {
    std::vector<Epoch> epochs;
    LogStats stats;
};

/// Incremental decoder: feed lines in order, collect epochs. A new epoch starts whenever a GGA or RMC
/// sentence carries a time tag different from the current one; GSA, GSV and VTG attach to the current
/// epoch.
class LogParser {
public:
    explicit LogParser(bool require_checksum = true) : require_checksum_(require_checksum) {}

    /// Decode one line. Returns the split status (ok for sentences of unsupported types too).
    Status feed(std::string_view line);

    const LogStats& stats() const { return stats_; }
    const std::vector<Epoch>& epochs() const { return epochs_; }
    Log finish() &&;

private:
    Epoch& epoch_for(const std::optional<double>& time_of_day_s);
    Epoch& current();

    bool require_checksum_;
    LogStats stats_;
    std::vector<Epoch> epochs_;
};

/// Decode a whole stream or file.
Log parse_log(std::istream& in, bool require_checksum = true);
Log parse_log_file(const std::string& path, bool require_checksum = true);

}  // namespace sensor_fusion::nmea
