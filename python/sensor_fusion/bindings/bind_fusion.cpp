// Chapter 6 bindings: the planar pose EKF and the GNSS + odometry fusion loop.
#include <pybind11/stl.h>

#include <limits>

#include "bindings.hpp"
#include "sensor_fusion/fusion.hpp"

namespace sensor_fusion::python {

namespace {
using RowsXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
}

void bind_fusion(py::module_ m) {
    using namespace fusion;
    m.def("wrap_angle", &wrap_angle, py::arg("a"));
    m.def("compose", &compose, py::arg("pose"), py::arg("delta"), "pose (+) body-frame increment, eq. (6.1)");
    m.def("between", &between, py::arg("a"), py::arg("b"), "body-frame increment from pose a to pose b");
    m.def("chi2_2dof", &chi2_2dof, py::arg("probability"));

    py::class_<OdometryNoise>(m, "OdometryNoise", "odometry increment covariance, eq. (6.3)")
        .def(py::init([](double kt, double kr, double krm, double ft, double fr) {
                 OdometryNoise n;
                 n.k_trans = kt, n.k_rot = kr, n.k_rot_per_m = krm, n.floor_trans = ft, n.floor_rot = fr;
                 return n;
             }),
             py::arg("k_trans") = 0.02, py::arg("k_rot") = 0.05, py::arg("k_rot_per_m") = 0.002,
             py::arg("floor_trans") = 1e-3, py::arg("floor_rot") = 1e-4)
        .def_readwrite("k_trans", &OdometryNoise::k_trans)
        .def_readwrite("k_rot", &OdometryNoise::k_rot)
        .def_readwrite("k_rot_per_m", &OdometryNoise::k_rot_per_m)
        .def_readwrite("floor_trans", &OdometryNoise::floor_trans)
        .def_readwrite("floor_rot", &OdometryNoise::floor_rot)
        .def("covariance", &OdometryNoise::covariance, py::arg("delta"));

    py::class_<PoseEkf>(m, "PoseEkf", "EKF on (x, y, yaw), eqs. (6.1)-(6.6)")
        .def(py::init<const Eigen::Vector3d&, const Eigen::Matrix3d&>(), py::arg("x0"), py::arg("P0"))
        .def("predict", &PoseEkf::predict, py::arg("delta"), py::arg("Q"))
        .def(
            "update_position",
            [](PoseEkf& e, const Eigen::Vector2d& z, const Eigen::Matrix2d& R, double gate) {
                const auto u = e.update_position(z, R, gate);
                return py::make_tuple(u.nis, u.accepted);
            },
            py::arg("z"), py::arg("R"), py::arg("gate") = std::numeric_limits<double>::infinity(),
            "returns (nis, accepted)")
        .def_property_readonly("x", &PoseEkf::state)
        .def_property_readonly("P", &PoseEkf::covariance);

    m.def(
        "fuse",
        [](const std::vector<double>& odom_times, const RowsXd& odom_poses, const std::vector<double>& gnss_times,
           const RowsXd& gnss_positions, const RowsXd& gnss_variances, const Eigen::Vector3d& x0,
           const Eigen::Matrix3d& P0, const OdometryNoise& noise, double gate_probability, double gnss_scale) {
            if (odom_poses.cols() != 3 || gnss_positions.cols() != 2 || gnss_variances.cols() != 2)
                throw std::invalid_argument("odom_poses (n, 3), gnss_positions (m, 2), gnss_variances (m, 2)");
            FusionInput in;
            in.odom_times = odom_times;
            in.gnss_times = gnss_times;
            for (Eigen::Index i = 0; i < odom_poses.rows(); ++i) in.odom_poses.emplace_back(odom_poses.row(i).transpose());
            for (Eigen::Index i = 0; i < gnss_positions.rows(); ++i) {
                in.gnss_positions.emplace_back(gnss_positions.row(i).transpose());
                in.gnss_covariances.emplace_back(Eigen::Vector2d(gnss_variances.row(i).transpose()).asDiagonal());
            }
            FusionOptions o;
            o.odometry = noise;
            o.gate_probability = gate_probability;
            o.gnss_scale = gnss_scale;
            const auto r = fuse(in, x0, P0, o);
            RowsXd X(r.states.size(), 3), S(r.states.size(), 3);
            for (std::size_t k = 0; k < r.states.size(); ++k) {
                X.row(Eigen::Index(k)) = r.states[k].transpose();
                S.row(Eigen::Index(k)) = r.covariances[k].diagonal().cwiseSqrt().transpose();
            }
            py::dict d;
            d["t"] = r.times;
            d["x"] = X;
            d["sigma"] = S;
            d["nis"] = r.nis;
            d["accepted"] = r.accepted;
            return d;
        },
        py::arg("odom_times"), py::arg("odom_poses"), py::arg("gnss_times"), py::arg("gnss_positions"),
        py::arg("gnss_variances"), py::arg("x0"), py::arg("P0"), py::arg("noise") = OdometryNoise{},
        py::arg("gate_probability") = 0.99, py::arg("gnss_scale") = 1.0,
        "Algorithm 6.1; GNSS covariances are diagonal (var_e, var_n). Returns {t, x, sigma, nis, accepted}");
}

}  // namespace sensor_fusion::python
