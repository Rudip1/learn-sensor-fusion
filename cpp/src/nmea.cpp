#include "sensor_fusion/nmea.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

namespace sensor_fusion::nmea {

namespace {

std::string_view trim(std::string_view s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.remove_suffix(1);
    return s;
}

int hex_value(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    return -1;
}

/// Strict decimal parse: the whole field must be a number.
std::optional<double> to_double(std::string_view s) {
    s = trim(s);
    if (s.empty()) return std::nullopt;
    const std::string buf(s);
    char* end = nullptr;
    const double v = std::strtod(buf.c_str(), &end);
    if (end != buf.c_str() + buf.size() || !std::isfinite(v)) return std::nullopt;
    return v;
}

std::optional<int> to_int(std::string_view s) {
    s = trim(s);
    if (s.empty()) return std::nullopt;
    const std::string buf(s);
    char* end = nullptr;
    const long v = std::strtol(buf.c_str(), &end, 10);
    if (end != buf.c_str() + buf.size()) return std::nullopt;
    return static_cast<int>(v);
}

std::string_view field(const Sentence& s, std::size_t i) {
    return i < s.fields.size() ? std::string_view(s.fields[i]) : std::string_view();
}

char first_char(std::string_view s) { return s.empty() ? 0 : s.front(); }

}  // namespace

const char* to_string(Status status) {
    switch (status) {
        case Status::ok: return "ok";
        case Status::empty: return "empty";
        case Status::no_start_delimiter: return "no_start_delimiter";
        case Status::missing_checksum: return "missing_checksum";
        case Status::bad_checksum_format: return "bad_checksum_format";
        case Status::checksum_mismatch: return "checksum_mismatch";
        case Status::malformed: return "malformed";
    }
    return "unknown";
}

std::uint8_t checksum(std::string_view body) {
    std::uint8_t c = 0;
    for (const char ch : body) c ^= static_cast<std::uint8_t>(ch);  // eq. (1.6)
    return c;
}

SplitResult split(std::string_view line, bool require_checksum) {
    SplitResult out;
    line = trim(line);
    if (line.empty()) return out;  // Status::empty

    const auto start = line.find('$');
    if (start == std::string_view::npos) {
        out.status = Status::no_start_delimiter;
        return out;
    }
    line.remove_prefix(start + 1);

    std::string_view body = line;
    const auto star = line.rfind('*');
    if (star == std::string_view::npos) {
        if (require_checksum) {
            out.status = Status::missing_checksum;
            return out;
        }
    } else {
        body = line.substr(0, star);
        const std::string_view hex = line.substr(star + 1);
        if (hex.size() != 2 || hex_value(hex[0]) < 0 || hex_value(hex[1]) < 0) {
            out.status = Status::bad_checksum_format;
            return out;
        }
        const int transmitted = hex_value(hex[0]) * 16 + hex_value(hex[1]);
        if (transmitted != checksum(body)) {
            out.status = Status::checksum_mismatch;
            return out;
        }
    }

    // Split at commas; the first token is the address field "TTSSS".
    std::vector<std::string> tokens;
    std::size_t pos = 0;
    while (true) {
        const auto comma = body.find(',', pos);
        tokens.emplace_back(body.substr(pos, comma == std::string_view::npos ? body.npos : comma - pos));
        if (comma == std::string_view::npos) break;
        pos = comma + 1;
    }
    const std::string& address = tokens.front();
    if (address.size() < 5) {
        out.status = Status::malformed;
        return out;
    }
    // Proprietary sentences ("$P...") have a one-letter prefix; standard ones a two-letter talker.
    if (address[0] == 'P') {
        out.sentence.talker = "P";
        out.sentence.type = address.substr(1);
    } else {
        out.sentence.talker = address.substr(0, 2);
        out.sentence.type = address.substr(2);
    }
    out.sentence.fields.assign(tokens.begin() + 1, tokens.end());
    out.status = Status::ok;
    return out;
}

std::optional<double> parse_angle(std::string_view value, std::string_view hemisphere) {
    const auto raw = to_double(value);
    const char h = first_char(trim(hemisphere));
    if (!raw || *raw < 0.0 || (h != 'N' && h != 'S' && h != 'E' && h != 'W')) return std::nullopt;
    // eq. (1.7): the integer hundreds are whole degrees, the remainder is minutes.
    const double degrees = std::floor(*raw / 100.0);
    const double minutes = *raw - 100.0 * degrees;
    if (minutes >= 60.0) return std::nullopt;
    const double angle = degrees + minutes / 60.0;
    return (h == 'S' || h == 'W') ? -angle : angle;
}

std::optional<double> parse_time_of_day(std::string_view value) {
    value = trim(value);
    if (value.size() < 6) return std::nullopt;
    const auto hh = to_int(value.substr(0, 2));
    const auto mm = to_int(value.substr(2, 2));
    const auto ss = to_double(value.substr(4));
    if (!hh || !mm || !ss || *hh > 23 || *mm > 59 || *ss >= 61.0) return std::nullopt;
    return 3600.0 * *hh + 60.0 * *mm + *ss;
}

std::optional<Date> parse_date(std::string_view value) {
    value = trim(value);
    if (value.size() != 6) return std::nullopt;
    const auto dd = to_int(value.substr(0, 2));
    const auto mo = to_int(value.substr(2, 2));
    const auto yy = to_int(value.substr(4, 2));
    if (!dd || !mo || !yy || *dd < 1 || *dd > 31 || *mo < 1 || *mo > 12) return std::nullopt;
    return Date{2000 + *yy, *mo, *dd};
}

std::optional<Gga> parse_gga(const Sentence& s) {
    // hhmmss.ss, lat, N/S, lon, E/W, quality, numSV, HDOP, alt, M, sep, M, age, station
    if (s.type != "GGA" || s.fields.size() < 9) return std::nullopt;
    Gga g;
    g.time_of_day_s = parse_time_of_day(field(s, 0));
    g.latitude_deg = parse_angle(field(s, 1), field(s, 2));
    g.longitude_deg = parse_angle(field(s, 3), field(s, 4));
    g.fix_quality = to_int(field(s, 5)).value_or(0);
    g.satellites_used = to_int(field(s, 6)).value_or(0);
    g.hdop = to_double(field(s, 7));
    g.altitude_msl_m = to_double(field(s, 8));
    g.geoid_separation_m = to_double(field(s, 10));
    g.dgps_age_s = to_double(field(s, 12));
    g.dgps_station = to_int(field(s, 13));
    return g;
}

std::optional<Rmc> parse_rmc(const Sentence& s) {
    // hhmmss.ss, status, lat, N/S, lon, E/W, SOG[kn], COG, ddmmyy, magvar, E/W, mode
    if (s.type != "RMC" || s.fields.size() < 9) return std::nullopt;
    Rmc r;
    r.time_of_day_s = parse_time_of_day(field(s, 0));
    r.valid = first_char(field(s, 1)) == 'A';
    r.latitude_deg = parse_angle(field(s, 2), field(s, 3));
    r.longitude_deg = parse_angle(field(s, 4), field(s, 5));
    if (const auto kn = to_double(field(s, 6))) r.speed_mps = knots_to_mps(*kn);
    r.course_deg = to_double(field(s, 7));
    r.date = parse_date(field(s, 8));
    r.mode = first_char(field(s, 11));
    return r;
}

std::optional<Gsa> parse_gsa(const Sentence& s) {
    // mode, fix type, 12 PRN slots, PDOP, HDOP, VDOP  (NMEA 4.1 appends a system ID)
    if (s.type != "GSA" || s.fields.size() < 17) return std::nullopt;
    Gsa g;
    g.selection_mode = first_char(field(s, 0));
    g.fix_type = to_int(field(s, 1)).value_or(1);
    for (std::size_t i = 2; i < 14; ++i)
        if (const auto prn = to_int(field(s, i))) g.prns.push_back(*prn);
    g.pdop = to_double(field(s, 14));
    g.hdop = to_double(field(s, 15));
    g.vdop = to_double(field(s, 16));
    return g;
}

std::optional<Gsv> parse_gsv(const Sentence& s) {
    // total, number, in view, then groups of (PRN, elevation, azimuth, SNR)
    if (s.type != "GSV" || s.fields.size() < 3) return std::nullopt;
    Gsv g;
    const auto total = to_int(field(s, 0));
    const auto number = to_int(field(s, 1));
    if (!total || !number) return std::nullopt;
    g.total_messages = *total;
    g.message_number = *number;
    g.satellites_in_view = to_int(field(s, 2)).value_or(0);
    // Only complete groups of four; NMEA 4.10 appends a lone signal-ID field that must not be read as a PRN.
    for (std::size_t i = 3; i + 3 < s.fields.size(); i += 4) {
        const auto prn = to_int(field(s, i));
        if (!prn) continue;
        SatelliteInView sat;
        sat.prn = *prn;
        sat.elevation_deg = to_double(field(s, i + 1));
        sat.azimuth_deg = to_double(field(s, i + 2));
        sat.snr_dbhz = to_double(field(s, i + 3));
        g.satellites.push_back(sat);
    }
    return g;
}

std::optional<Vtg> parse_vtg(const Sentence& s) {
    // COG true, T, COG magnetic, M, SOG kn, N, SOG km/h, K, mode
    if (s.type != "VTG" || s.fields.size() < 8) return std::nullopt;
    Vtg v;
    v.course_true_deg = to_double(field(s, 0));
    v.course_magnetic_deg = to_double(field(s, 2));
    if (const auto kmh = to_double(field(s, 6)))
        v.speed_mps = *kmh / 3.6;
    else if (const auto kn = to_double(field(s, 4)))
        v.speed_mps = knots_to_mps(*kn);
    v.mode = first_char(field(s, 8));
    return v;
}

Epoch& LogParser::current() {
    if (epochs_.empty()) epochs_.emplace_back();
    return epochs_.back();
}

Epoch& LogParser::epoch_for(const std::optional<double>& t) {
    if (epochs_.empty()) {
        epochs_.emplace_back();
        epochs_.back().time_of_day_s = t;
        return epochs_.back();
    }
    Epoch& last = epochs_.back();
    if (!t) return last;
    if (!last.time_of_day_s) {
        last.time_of_day_s = t;  // the epoch was opened by an untimed sentence
        return last;
    }
    if (std::abs(*last.time_of_day_s - *t) < 1e-6) return last;
    epochs_.emplace_back();
    epochs_.back().time_of_day_s = t;
    return epochs_.back();
}

Status LogParser::feed(std::string_view line) {
    ++stats_.lines;
    const SplitResult r = split(line, require_checksum_);
    if (r.status == Status::empty) return r.status;
    if (r.status == Status::checksum_mismatch || r.status == Status::bad_checksum_format ||
        r.status == Status::missing_checksum) {
        ++stats_.checksum_errors;
        return r.status;
    }
    if (r.status != Status::ok) {
        ++stats_.other_errors;
        return r.status;
    }
    ++stats_.sentences;
    const Sentence& s = r.sentence;
    ++stats_.counts_by_type[s.type];

    if (s.type == "GGA") {
        const auto g = parse_gga(s);
        if (!g) return ++stats_.other_errors, Status::malformed;
        Epoch& e = epoch_for(g->time_of_day_s);
        e.fix_quality = g->fix_quality;
        e.satellites_used = g->satellites_used;
        e.latitude_deg = g->latitude_deg;
        e.longitude_deg = g->longitude_deg;
        e.altitude_msl_m = g->altitude_msl_m;
        e.geoid_separation_m = g->geoid_separation_m;
        e.hdop = g->hdop;
    } else if (s.type == "RMC") {
        const auto m = parse_rmc(s);
        if (!m) return ++stats_.other_errors, Status::malformed;
        Epoch& e = epoch_for(m->time_of_day_s);
        e.rmc_valid = m->valid;
        e.date = m->date;
        e.speed_mps = m->speed_mps;
        e.course_deg = m->course_deg;
        if (!e.latitude_deg) e.latitude_deg = m->latitude_deg;
        if (!e.longitude_deg) e.longitude_deg = m->longitude_deg;
    } else if (s.type == "GSA") {
        const auto g = parse_gsa(s);
        if (!g) return ++stats_.other_errors, Status::malformed;
        Epoch& e = current();
        e.fix_type = g->fix_type;
        e.used_prns = g->prns;
        e.pdop = g->pdop;
        e.vdop = g->vdop;
        if (!e.hdop) e.hdop = g->hdop;
    } else if (s.type == "GSV") {
        const auto g = parse_gsv(s);
        if (!g) return ++stats_.other_errors, Status::malformed;
        Epoch& e = current();
        if (g->message_number == 1) e.satellites.clear();
        e.satellites.insert(e.satellites.end(), g->satellites.begin(), g->satellites.end());
    } else if (s.type == "VTG") {
        const auto v = parse_vtg(s);
        if (!v) return ++stats_.other_errors, Status::malformed;
        Epoch& e = current();
        if (!e.speed_mps) e.speed_mps = v->speed_mps;
        if (!e.course_deg) e.course_deg = v->course_true_deg;
    }
    return Status::ok;
}

Log LogParser::finish() && { return Log{std::move(epochs_), std::move(stats_)}; }

Log parse_log(std::istream& in, bool require_checksum) {
    LogParser parser(require_checksum);
    std::string line;
    while (std::getline(in, line)) parser.feed(line);
    return std::move(parser).finish();
}

Log parse_log_file(const std::string& path, bool require_checksum) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open NMEA log: " + path);
    return parse_log(in, require_checksum);
}

}  // namespace sensor_fusion::nmea
