#include <pybind11/stl.h>

#include <cstdio>

#include "bindings.hpp"
#include "sensor_fusion/geodesy.hpp"
#include "sensor_fusion/gnss.hpp"
#include "sensor_fusion/gnss_quality.hpp"

namespace sensor_fusion::python {

void bind_gnss(py::module_ m) {
    using namespace sensor_fusion::gnss;

    m.attr("SPEED_OF_LIGHT") = kSpeedOfLight;

    py::class_<PositionSolution>(m, "PositionSolution")
        .def_readonly("position", &PositionSolution::position, "receiver position, ECEF [m]")
        .def_readonly("clock_bias_m", &PositionSolution::clock_bias_m)
        .def_readonly("cofactor", &PositionSolution::cofactor, "(G^T W G)^-1, eq. (1.4)")
        .def_readonly("residuals", &PositionSolution::residuals)
        .def_readonly("iterations", &PositionSolution::iterations)
        .def_readonly("converged", &PositionSolution::converged);

    m.def(
        "geometry_matrix",
        [](const Eigen::Vector3d& receiver, const RowsX3& satellites) {
            return Eigen::MatrixX4d(geometry_matrix(receiver, satellites.transpose()));
        },
        py::arg("receiver"), py::arg("satellites"), "rows [-u_i^T, 1], eq. (1.2); satellites as (n, 3)");
    m.def(
        "geometric_ranges",
        [](const Eigen::Vector3d& receiver, const RowsX3& satellites) {
            return Eigen::VectorXd(geometric_ranges(receiver, satellites.transpose()));
        },
        py::arg("receiver"), py::arg("satellites"));
    m.def(
        "solve_position",
        [](const RowsX3& satellites, const Eigen::VectorXd& pseudoranges, std::optional<Eigen::VectorXd> weights,
           std::optional<Eigen::Vector3d> initial_position, double initial_clock_bias_m, int max_iterations,
           double tolerance_m) {
            SolverOptions opt;
            if (initial_position) opt.initial_position = *initial_position;
            opt.initial_clock_bias_m = initial_clock_bias_m;
            opt.max_iterations = max_iterations;
            opt.tolerance_m = tolerance_m;
            const Eigen::Matrix3Xd s = satellites.transpose();
            return solve_position(s, pseudoranges, weights.value_or(Eigen::VectorXd()), opt);
        },
        py::arg("satellites"), py::arg("pseudoranges"), py::arg("weights") = py::none(),
        py::arg("initial_position") = py::none(), py::arg("initial_clock_bias_m") = 0.0,
        py::arg("max_iterations") = 20, py::arg("tolerance_m") = 1e-4,
        "Gauss-Newton solution of rho_i = |s_i - r| + b, eqs. (1.1)-(1.3). satellites: (n, 3) ECEF [m]");
    m.def(
        "place_satellites",
        [](const Eigen::Vector3d& receiver, const Eigen::VectorXd& azimuth, const Eigen::VectorXd& elevation,
           double orbit_radius) {
            return RowsX3(place_satellites(receiver, azimuth, elevation, orbit_radius).transpose());
        },
        py::arg("receiver"), py::arg("azimuth"), py::arg("elevation"), py::arg("orbit_radius") = 26'560'000.0,
        "satellites (n, 3) at the given azimuth/elevation [rad], eq. (1.5)");

    // ---- chapter 3: quality -------------------------------------------------------------------------
    using geo::deg2rad;
    auto rad = [](const Eigen::VectorXd& deg) { return Eigen::VectorXd(deg * (3.14159265358979323846 / 180.0)); };

    py::class_<Dop>(m, "Dop")
        .def_readonly("gdop", &Dop::gdop)
        .def_readonly("pdop", &Dop::pdop)
        .def_readonly("hdop", &Dop::hdop)
        .def_readonly("vdop", &Dop::vdop)
        .def_readonly("tdop", &Dop::tdop)
        .def_readonly("cofactor_enu", &Dop::cofactor_enu)
        .def_readonly("valid", &Dop::valid)
        .def("__repr__", [](const Dop& d) {
            char buf[160];
            std::snprintf(buf, sizeof buf, "Dop(gdop=%.3f, pdop=%.3f, hdop=%.3f, vdop=%.3f, tdop=%.3f, valid=%s)",
                          d.gdop, d.pdop, d.hdop, d.vdop, d.tdop, d.valid ? "True" : "False");
            return std::string(buf);
        });

    m.def(
        "geometry_matrix_enu",
        [rad](const Eigen::VectorXd& az_deg, const Eigen::VectorXd& el_deg) {
            return Eigen::MatrixX4d(geometry_matrix_enu(rad(az_deg), rad(el_deg)));
        },
        py::arg("azimuth_deg"), py::arg("elevation_deg"), "rows [-u^T, 1] in ENU, eq. (3.1)");
    m.def(
        "compute_dop",
        [rad](const Eigen::VectorXd& az_deg, const Eigen::VectorXd& el_deg, std::optional<Eigen::VectorXd> w) {
            return compute_dop(rad(az_deg), rad(el_deg), w.value_or(Eigen::VectorXd()));
        },
        py::arg("azimuth_deg"), py::arg("elevation_deg"), py::arg("weights") = py::none(), "eqs. (3.2)-(3.3)");
    m.def(
        "dop_from_ecef_cofactor",
        [](const Eigen::Matrix4d& Q, double lat_deg, double lon_deg) {
            return dop_from_ecef_cofactor(Q, deg2rad(lat_deg), deg2rad(lon_deg));
        },
        py::arg("cofactor_ecef"), py::arg("lat_deg"), py::arg("lon_deg"), "eq. (3.4)");
    m.def(
        "elevation_weights", [rad](const Eigen::VectorXd& el_deg) { return elevation_weights(rad(el_deg)); },
        py::arg("elevation_deg"), "w = sin^2(el)");

    py::class_<ErrorEllipse>(m, "ErrorEllipse")
        .def_readonly("semi_major", &ErrorEllipse::semi_major)
        .def_readonly("semi_minor", &ErrorEllipse::semi_minor)
        .def_readonly("angle", &ErrorEllipse::angle, "major axis direction, counter-clockwise from east [rad]")
        .def("__repr__", [](const ErrorEllipse& e) {
            char buf[120];
            std::snprintf(buf, sizeof buf, "ErrorEllipse(semi_major=%.4f, semi_minor=%.4f, angle_deg=%.2f)",
                          e.semi_major, e.semi_minor, e.angle * 57.29577951308232);
            return std::string(buf);
        });
    m.def("error_ellipse", &error_ellipse, py::arg("cov"), py::arg("probability") = 0.3934693402873666,
          "confidence ellipse of a 2-D Gaussian, eq. (3.7)");

    py::class_<PrecisionStats>(m, "PrecisionStats")
        .def_readonly("mean", &PrecisionStats::mean)
        .def_readonly("covariance", &PrecisionStats::covariance)
        .def_readonly("sigma_east", &PrecisionStats::sigma_east)
        .def_readonly("sigma_north", &PrecisionStats::sigma_north)
        .def_readonly("drms", &PrecisionStats::drms)
        .def_readonly("two_drms", &PrecisionStats::two_drms)
        .def_readonly("cep50", &PrecisionStats::cep50)
        .def_readonly("r95", &PrecisionStats::r95)
        .def_readonly("count", &PrecisionStats::count);
    m.def(
        "precision_stats", [](const RowsX2& en) { return precision_stats(en.transpose()); }, py::arg("en"),
        "precision of (n, 2) east/north positions, eq. (3.6)");

    m.def("gga_quality_name", &gga_quality_name, py::arg("quality"));
    m.def("gsa_fix_type_name", &gsa_fix_type_name, py::arg("fix_type"));
    m.def("navsat_status_name", &navsat_status_name, py::arg("status"));

    py::class_<MotionParams>(m, "MotionParams")
        .def(py::init<>())
        .def_readwrite("window_s", &MotionParams::window_s)
        .def_readwrite("sigma_displacement_m", &MotionParams::sigma_displacement_m)
        .def_readwrite("chi2_threshold", &MotionParams::chi2_threshold)
        .def_readwrite("speed_threshold_mps", &MotionParams::speed_threshold_mps)
        .def_readwrite("min_run", &MotionParams::min_run);
    py::class_<MotionResult>(m, "MotionResult")
        .def_readonly("moving", &MotionResult::moving)
        .def_readonly("test_statistic", &MotionResult::test_statistic)
        .def_readonly("displacement_speed", &MotionResult::displacement_speed);
    m.def(
        "classify_motion",
        [](const Eigen::VectorXd& t, const RowsX2& en, std::optional<Eigen::VectorXd> speed, const MotionParams& p) {
            return classify_motion(t, en.transpose(), speed.value_or(Eigen::VectorXd()), p);
        },
        py::arg("t"), py::arg("en"), py::arg("speed") = py::none(), py::arg("params") = MotionParams(),
        "Algorithm 3.1: t (n,), en (n, 2), speed (n,) or None");
}

}  // namespace sensor_fusion::python
