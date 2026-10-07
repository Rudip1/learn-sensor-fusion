// Checks of the NMEA decoder against the hand-worked examples in 1_theory/01_gnss_fundamentals.md and
// against the recorded logs in data/gnss/.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <sstream>
#include <string>

#include "sensor_fusion/nmea.hpp"

using Catch::Approx;
namespace nmea = sensor_fusion::nmea;

namespace {
const char* kGga = "$GPGGA,091509.000,4728.4639,N,01903.4171,E,1,7,1.55,105.5,M,41.1,M,,*52";
const char* kGgaNoFix = "$GPGGA,090550.000,,,,,0,2,,,M,,M,,*43";
const char* kRmc = "$GPRMC,091509.000,A,4728.4639,N,01903.4171,E,2.13,143.84,041125,,,A*6A";
const char* kGsa = "$GPGSA,A,3,05,09,21,30,07,11,20,,,,,,1.79,1.55,0.89*04";
const char* kGsaNoFix = "$GPGSA,A,1,,,,,,,,,,,,,,,*1E";
const char* kGsvLast = "$GPGSV,4,4,14,04,02,097,,27,02,055,*73";
const char* kVtg = "$GPVTG,143.84,T,,M,2.13,N,3.94,K,A*39";
}  // namespace

TEST_CASE("checksum: worked examples of eq. (1.6)", "[nmea][ch1]") {
    CHECK(nmea::checksum("GPGSA,A,1,,,,,,,,,,,,,,,") == 0x1E);
    CHECK(nmea::checksum("GPGGA,091509.000,4728.4639,N,01903.4171,E,1,7,1.55,105.5,M,41.1,M,,") == 0x52);
    CHECK(nmea::checksum("") == 0x00);
    CHECK(nmea::checksum("AA") == 0x00);  // equal bytes cancel
}

TEST_CASE("split: address, fields and checksum verification", "[nmea][ch1]") {
    const auto r = nmea::split(kGga);
    REQUIRE(r.status == nmea::Status::ok);
    CHECK(r.sentence.talker == "GP");
    CHECK(r.sentence.type == "GGA");
    REQUIRE(r.sentence.fields.size() == 14);
    CHECK(r.sentence.fields[0] == "091509.000");
    CHECK(r.sentence.fields[12].empty());

    std::string corrupted = kGga;
    corrupted[20] = '5';  // one digit of the latitude
    CHECK(nmea::split(corrupted).status == nmea::Status::checksum_mismatch);
    CHECK(nmea::split("$GPGGA,1,2,3").status == nmea::Status::missing_checksum);
    CHECK(nmea::split("$GPGGA,1,2,3", false).status == nmea::Status::ok);
    CHECK(nmea::split("$GPGGA,1*5").status == nmea::Status::bad_checksum_format);
    CHECK(nmea::split("$GPGGA,1*ZZ").status == nmea::Status::bad_checksum_format);
    CHECK(nmea::split("GPGGA,1*00").status == nmea::Status::no_start_delimiter);
    CHECK(nmea::split("   ").status == nmea::Status::empty);
    CHECK(nmea::split("$GP*17").status == nmea::Status::malformed);
    // leading noise and a trailing CR are tolerated
    CHECK(nmea::split(std::string("}") + kGsa + "\r").status == nmea::Status::ok);
    // lower-case hex digits are accepted
    CHECK(nmea::split("$GPGSA,A,1,,,,,,,,,,,,,,,*1e").status == nmea::Status::ok);
}

TEST_CASE("parse_angle: ddmm.mmmm to signed degrees, eq. (1.7)", "[nmea][ch1]") {
    CHECK(*nmea::parse_angle("4728.4639", "N") == Approx(47.474398333333).epsilon(1e-13));
    CHECK(*nmea::parse_angle("01903.4171", "E") == Approx(19.056951666667).epsilon(1e-13));
    CHECK(*nmea::parse_angle("4728.4639", "S") == Approx(-47.474398333333).epsilon(1e-13));
    CHECK(*nmea::parse_angle("12000.0000", "W") == Approx(-120.0));
    CHECK_FALSE(nmea::parse_angle("", "N"));
    CHECK_FALSE(nmea::parse_angle("4728.4639", ""));
    CHECK_FALSE(nmea::parse_angle("4775.0000", "N"));  // 75 minutes is not an angle
    CHECK_FALSE(nmea::parse_angle("47x8.4639", "N"));
}

TEST_CASE("time and date fields", "[nmea][ch1]") {
    CHECK(*nmea::parse_time_of_day("091509.000") == Approx(33309.0));
    CHECK(*nmea::parse_time_of_day("090554.202") == Approx(32754.202));
    CHECK_FALSE(nmea::parse_time_of_day("2515"));
    const auto d = nmea::parse_date("041125");
    REQUIRE(d);
    CHECK(d->year == 2025);
    CHECK(d->month == 11);
    CHECK(d->day == 4);
    CHECK_FALSE(nmea::parse_date("321325"));
}

TEST_CASE("GGA: worked example and a sentence without fix", "[nmea][ch1]") {
    const auto g = nmea::parse_gga(nmea::split(kGga).sentence);
    REQUIRE(g);
    CHECK(*g->time_of_day_s == Approx(33309.0));
    CHECK(*g->latitude_deg == Approx(47.474398333333));
    CHECK(*g->longitude_deg == Approx(19.056951666667));
    CHECK(g->fix_quality == 1);
    CHECK(g->satellites_used == 7);
    CHECK(*g->hdop == Approx(1.55));
    CHECK(*g->altitude_msl_m == Approx(105.5));
    CHECK(*g->geoid_separation_m == Approx(41.1));
    CHECK(*g->altitude_msl_m + *g->geoid_separation_m == Approx(146.6));  // eq. (1.8)
    CHECK_FALSE(g->dgps_age_s);

    const auto n = nmea::parse_gga(nmea::split(kGgaNoFix).sentence);
    REQUIRE(n);
    CHECK(n->fix_quality == 0);
    CHECK(n->satellites_used == 2);
    CHECK_FALSE(n->latitude_deg);
    CHECK_FALSE(n->hdop);
    CHECK_FALSE(nmea::parse_gga(nmea::split(kRmc).sentence));  // wrong type
}

TEST_CASE("RMC, GSA, GSV and VTG", "[nmea][ch1]") {
    const auto r = nmea::parse_rmc(nmea::split(kRmc).sentence);
    REQUIRE(r);
    CHECK(r->valid);
    CHECK(*r->speed_mps == Approx(2.13 * 1852.0 / 3600.0));  // eq. (1.9)
    CHECK(*r->course_deg == Approx(143.84));
    CHECK(r->date->day == 4);
    CHECK(r->mode == 'A');

    const auto a = nmea::parse_gsa(nmea::split(kGsa).sentence);
    REQUIRE(a);
    CHECK(a->fix_type == 3);
    CHECK(a->prns == std::vector<int>{5, 9, 21, 30, 7, 11, 20});
    CHECK(*a->pdop == Approx(1.79));
    CHECK(*a->hdop == Approx(1.55));
    CHECK(*a->vdop == Approx(0.89));
    const auto a0 = nmea::parse_gsa(nmea::split(kGsaNoFix).sentence);
    REQUIRE(a0);
    CHECK(a0->fix_type == 1);
    CHECK(a0->prns.empty());
    CHECK_FALSE(a0->pdop);

    const auto v = nmea::parse_gsv(nmea::split(kGsvLast).sentence);
    REQUIRE(v);
    CHECK(v->total_messages == 4);
    CHECK(v->message_number == 4);
    CHECK(v->satellites_in_view == 14);
    REQUIRE(v->satellites.size() == 2);
    CHECK(v->satellites[0].prn == 4);
    CHECK(*v->satellites[0].elevation_deg == Approx(2.0));
    CHECK(*v->satellites[0].azimuth_deg == Approx(97.0));
    CHECK_FALSE(v->satellites[0].snr_dbhz);

    const auto t = nmea::parse_vtg(nmea::split(kVtg).sentence);
    REQUIRE(t);
    CHECK(*t->course_true_deg == Approx(143.84));
    CHECK_FALSE(t->course_magnetic_deg);
    CHECK(*t->speed_mps == Approx(3.94 / 3.6));
}

TEST_CASE("LogParser groups sentences into epochs", "[nmea][ch1]") {
    std::istringstream in(std::string(kGga) + "\n" + kGsa + "\n" +
                          "$GPGSV,2,1,05,30,66,218,33,09,42,091,24,05,31,311,19,16,09,032,*75\n"
                          "$GPGSV,2,2,05,21,,,23*7E\n" +
                          kRmc + "\n" + kVtg + "\n" + "$GPGGA,091510.000,,,,,0,0,,,M,,M,,*44\n" + "garbage\n");
    const nmea::Log log = nmea::parse_log(in);
    REQUIRE(log.epochs.size() == 2);
    const auto& e = log.epochs[0];
    CHECK(*e.time_of_day_s == Approx(33309.0));
    CHECK(e.fix_type == 3);
    CHECK(e.used_prns.size() == 7);
    CHECK(e.satellites.size() == 5);
    CHECK(*e.satellites[4].snr_dbhz == Approx(23.0));
    CHECK(e.rmc_valid);
    CHECK(*e.speed_mps == Approx(2.13 * 1852.0 / 3600.0));  // RMC wins over VTG
    CHECK(e.date->month == 11);
    CHECK(log.stats.sentences == 7);
    CHECK(log.stats.other_errors == 1);
    CHECK(log.stats.counts_by_type.at("GSV") == 2);
}

TEST_CASE("recorded logs decode without checksum errors", "[nmea][ch1][data]") {
    for (const char* name : {"nmea_log1.nmea", "nmea_log2.nmea"}) {
        const auto log = nmea::parse_log_file(std::string(SF_DATA_DIR) + "/gnss/" + name);
        INFO(name);
        CHECK(log.stats.checksum_errors == 0);
        CHECK(log.stats.other_errors == 0);
        CHECK(log.stats.sentences == log.stats.lines);
        CHECK(log.epochs.size() == static_cast<std::size_t>(log.stats.counts_by_type.at("GGA")));
    }
    const auto log2 = nmea::parse_log_file(std::string(SF_DATA_DIR) + "/gnss/nmea_log2.nmea");
    // the four-part GSV group at 09:15:13 lists all fourteen satellites in view
    bool found = false;
    for (const auto& e : log2.epochs) {
        if (e.time_of_day_s && std::abs(*e.time_of_day_s - (9 * 3600 + 15 * 60 + 13)) < 1e-6) {
            found = true;
            CHECK(e.satellites.size() == 14);
        }
    }
    CHECK(found);
}
