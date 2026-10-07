#include <pybind11/stl.h>

#include <cstdio>

#include "bindings.hpp"
#include "sensor_fusion/nmea.hpp"

namespace sensor_fusion::python {

void bind_nmea(py::module_ m) {
    using namespace sensor_fusion::nmea;

    py::enum_<Status>(m, "Status")
        .value("ok", Status::ok)
        .value("empty", Status::empty)
        .value("no_start_delimiter", Status::no_start_delimiter)
        .value("missing_checksum", Status::missing_checksum)
        .value("bad_checksum_format", Status::bad_checksum_format)
        .value("checksum_mismatch", Status::checksum_mismatch)
        .value("malformed", Status::malformed);

    py::class_<Sentence>(m, "Sentence")
        .def_readonly("talker", &Sentence::talker)
        .def_readonly("type", &Sentence::type)
        .def_readonly("fields", &Sentence::fields)
        .def("__repr__", [](const Sentence& s) {
            return "<Sentence " + s.talker + s.type + " with " + std::to_string(s.fields.size()) + " fields>";
        });

    py::class_<SplitResult>(m, "SplitResult")
        .def_readonly("status", &SplitResult::status)
        .def_readonly("sentence", &SplitResult::sentence);

    py::class_<Date>(m, "Date")
        .def_readonly("year", &Date::year)
        .def_readonly("month", &Date::month)
        .def_readonly("day", &Date::day)
        .def("__repr__", [](const Date& d) {
            char buf[32];
            std::snprintf(buf, sizeof buf, "%04d-%02d-%02d", d.year, d.month, d.day);
            return std::string(buf);
        });

    m.def("checksum", &checksum, py::arg("body"), "XOR of the characters between '$' and '*', eq. (1.6)");
    m.def("split", &split, py::arg("line"), py::arg("require_checksum") = true);
    m.def("parse_angle", &parse_angle, py::arg("value"), py::arg("hemisphere"),
          "'ddmm.mmmm' + hemisphere -> signed decimal degrees, eq. (1.7)");
    m.def("parse_time_of_day", &parse_time_of_day, py::arg("value"));
    m.def("parse_date", &parse_date, py::arg("value"));
    m.def("knots_to_mps", &knots_to_mps, py::arg("knots"));

    py::class_<Gga>(m, "Gga")
        .def_readonly("time_of_day_s", &Gga::time_of_day_s)
        .def_readonly("latitude_deg", &Gga::latitude_deg)
        .def_readonly("longitude_deg", &Gga::longitude_deg)
        .def_readonly("fix_quality", &Gga::fix_quality)
        .def_readonly("satellites_used", &Gga::satellites_used)
        .def_readonly("hdop", &Gga::hdop)
        .def_readonly("altitude_msl_m", &Gga::altitude_msl_m)
        .def_readonly("geoid_separation_m", &Gga::geoid_separation_m)
        .def_readonly("dgps_age_s", &Gga::dgps_age_s)
        .def_readonly("dgps_station", &Gga::dgps_station);
    py::class_<Rmc>(m, "Rmc")
        .def_readonly("time_of_day_s", &Rmc::time_of_day_s)
        .def_readonly("valid", &Rmc::valid)
        .def_readonly("latitude_deg", &Rmc::latitude_deg)
        .def_readonly("longitude_deg", &Rmc::longitude_deg)
        .def_readonly("speed_mps", &Rmc::speed_mps)
        .def_readonly("course_deg", &Rmc::course_deg)
        .def_readonly("date", &Rmc::date)
        .def_readonly("mode", &Rmc::mode);
    py::class_<Gsa>(m, "Gsa")
        .def_readonly("selection_mode", &Gsa::selection_mode)
        .def_readonly("fix_type", &Gsa::fix_type)
        .def_readonly("prns", &Gsa::prns)
        .def_readonly("pdop", &Gsa::pdop)
        .def_readonly("hdop", &Gsa::hdop)
        .def_readonly("vdop", &Gsa::vdop);
    py::class_<SatelliteInView>(m, "SatelliteInView")
        .def_readonly("prn", &SatelliteInView::prn)
        .def_readonly("elevation_deg", &SatelliteInView::elevation_deg)
        .def_readonly("azimuth_deg", &SatelliteInView::azimuth_deg)
        .def_readonly("snr_dbhz", &SatelliteInView::snr_dbhz);
    py::class_<Gsv>(m, "Gsv")
        .def_readonly("total_messages", &Gsv::total_messages)
        .def_readonly("message_number", &Gsv::message_number)
        .def_readonly("satellites_in_view", &Gsv::satellites_in_view)
        .def_readonly("satellites", &Gsv::satellites);
    py::class_<Vtg>(m, "Vtg")
        .def_readonly("course_true_deg", &Vtg::course_true_deg)
        .def_readonly("course_magnetic_deg", &Vtg::course_magnetic_deg)
        .def_readonly("speed_mps", &Vtg::speed_mps)
        .def_readonly("mode", &Vtg::mode);

    m.def("parse_gga", &parse_gga, py::arg("sentence"));
    m.def("parse_rmc", &parse_rmc, py::arg("sentence"));
    m.def("parse_gsa", &parse_gsa, py::arg("sentence"));
    m.def("parse_gsv", &parse_gsv, py::arg("sentence"));
    m.def("parse_vtg", &parse_vtg, py::arg("sentence"));

    py::class_<Epoch>(m, "Epoch")
        .def_readonly("time_of_day_s", &Epoch::time_of_day_s)
        .def_readonly("date", &Epoch::date)
        .def_readonly("fix_quality", &Epoch::fix_quality)
        .def_readonly("satellites_used", &Epoch::satellites_used)
        .def_readonly("latitude_deg", &Epoch::latitude_deg)
        .def_readonly("longitude_deg", &Epoch::longitude_deg)
        .def_readonly("altitude_msl_m", &Epoch::altitude_msl_m)
        .def_readonly("geoid_separation_m", &Epoch::geoid_separation_m)
        .def_readonly("hdop", &Epoch::hdop)
        .def_readonly("fix_type", &Epoch::fix_type)
        .def_readonly("used_prns", &Epoch::used_prns)
        .def_readonly("pdop", &Epoch::pdop)
        .def_readonly("vdop", &Epoch::vdop)
        .def_readonly("rmc_valid", &Epoch::rmc_valid)
        .def_readonly("speed_mps", &Epoch::speed_mps)
        .def_readonly("course_deg", &Epoch::course_deg)
        .def_readonly("satellites", &Epoch::satellites);

    py::class_<LogStats>(m, "LogStats")
        .def_readonly("lines", &LogStats::lines)
        .def_readonly("sentences", &LogStats::sentences)
        .def_readonly("checksum_errors", &LogStats::checksum_errors)
        .def_readonly("other_errors", &LogStats::other_errors)
        .def_readonly("counts_by_type", &LogStats::counts_by_type);

    py::class_<Log>(m, "Log").def_readonly("epochs", &Log::epochs).def_readonly("stats", &Log::stats);

    py::class_<LogParser>(m, "LogParser")
        .def(py::init<bool>(), py::arg("require_checksum") = true)
        .def("feed", &LogParser::feed, py::arg("line"))
        .def_property_readonly("stats", &LogParser::stats)
        .def_property_readonly("epochs", &LogParser::epochs);

    m.def("parse_log_file", &parse_log_file, py::arg("path"), py::arg("require_checksum") = true);
    m.def(
        "parse_log_lines",
        [](const std::vector<std::string>& lines, bool require_checksum) {
            LogParser p(require_checksum);
            for (const auto& l : lines) p.feed(l);
            return std::move(p).finish();
        },
        py::arg("lines"), py::arg("require_checksum") = true);
}

}  // namespace sensor_fusion::python
