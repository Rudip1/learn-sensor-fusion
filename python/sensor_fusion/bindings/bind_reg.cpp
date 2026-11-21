// Chapter 5 bindings: best_rigid, ICP and LiDAR odometry. Points are (n, 2) or (n, 3) arrays; transforms are
// homogeneous (3 x 3 or 4 x 4) matrices.
#include <pybind11/stl.h>

#include <optional>

#include "bindings.hpp"
#include "sensor_fusion/registration.hpp"

namespace sensor_fusion::python {

namespace {

using RowsXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template <int Dim>
cloud::Points<Dim> cols(const RowsXd& r) {
    if (r.cols() != Dim) throw std::invalid_argument("points must have shape (n, " + std::to_string(Dim) + ")");
    return r.transpose();
}

template <int Dim>
reg::Rigid<Dim> rigid(const Eigen::MatrixXd& M) {
    if (M.rows() != Dim + 1 || M.cols() != Dim + 1) throw std::invalid_argument("transform has the wrong size");
    reg::Rigid<Dim> T = reg::Rigid<Dim>::Identity();
    T.linear() = M.topLeftCorner(Dim, Dim);
    T.translation() = M.topRightCorner(Dim, 1);
    return T;
}

reg::IcpOptions options(const std::string& metric, int max_iterations, double max_correspondence, double trim,
                        double tolerance, int normal_k) {
    reg::IcpOptions o;
    if (metric == "point_to_point")
        o.metric = reg::Metric::PointToPoint;
    else if (metric == "point_to_plane")
        o.metric = reg::Metric::PointToPlane;
    else
        throw std::invalid_argument("metric must be 'point_to_point' or 'point_to_plane'");
    o.max_iterations = max_iterations;
    o.max_correspondence = max_correspondence;
    o.trim = trim;
    o.tolerance = tolerance;
    o.normal_k = normal_k;
    return o;
}

template <int Dim>
py::dict result(const reg::IcpResult<Dim>& r) {
    py::dict d;
    d["T"] = Eigen::MatrixXd(r.T.matrix());
    d["iterations"] = r.iterations;
    d["converged"] = r.converged;
    d["rmse"] = r.rmse;
    d["pairs"] = r.pairs;
    d["history"] = r.history;
    return d;
}

}  // namespace

void bind_reg(py::module_ m) {
    m.def(
        "best_rigid",
        [](const RowsXd& src, const RowsXd& dst) -> Eigen::MatrixXd {
            if (src.cols() == 2) return reg::best_rigid<2>(cols<2>(src), cols<2>(dst)).matrix();
            return reg::best_rigid<3>(cols<3>(src), cols<3>(dst)).matrix();
        },
        py::arg("src"), py::arg("dst"), "rigid T minimizing sum |dst - T src|^2 for matched rows, eqs. (5.2)-(5.3)");
    m.def(
        "icp",
        [](const RowsXd& source, const RowsXd& target, const std::optional<Eigen::MatrixXd>& initial,
           const std::string& metric, int max_iterations, double max_correspondence, double trim, double tolerance,
           int normal_k) -> py::dict {
            const auto o = options(metric, max_iterations, max_correspondence, trim, tolerance, normal_k);
            if (source.cols() == 2)
                return result<2>(reg::icp<2>(cols<2>(source), cols<2>(target),
                                             initial ? rigid<2>(*initial) : reg::Rigid<2>::Identity(), o));
            return result<3>(reg::icp<3>(cols<3>(source), cols<3>(target),
                                         initial ? rigid<3>(*initial) : reg::Rigid<3>::Identity(), o));
        },
        py::arg("source"), py::arg("target"), py::arg("initial") = py::none(), py::arg("metric") = "point_to_point",
        py::arg("max_iterations") = 50, py::arg("max_correspondence") = 1.0, py::arg("trim") = 1.0,
        py::arg("tolerance") = 1e-6, py::arg("normal_k") = 10,
        "ICP (Algorithm 5.1); returns {T, iterations, converged, rmse, pairs, history}; T maps source into target");
    m.def(
        "lidar_odometry",
        [](const std::vector<RowsXd>& scans, const std::string& metric, int max_iterations, double max_correspondence,
           double trim, double tolerance, int normal_k, double voxel, bool constant_velocity) -> py::dict {
            reg::OdometryOptions o;
            o.icp = options(metric, max_iterations, max_correspondence, trim, tolerance, normal_k);
            o.voxel = voxel;
            o.constant_velocity = constant_velocity;
            std::vector<Eigen::MatrixXd> poses;
            std::vector<double> rmse;
            std::vector<int> iterations;
            auto run = [&](auto dim) {
                constexpr int D = decltype(dim)::value;
                std::vector<cloud::Points<D>> s;
                for (const auto& r : scans) s.push_back(cols<D>(r));
                const auto res = reg::lidar_odometry<D>(s, o);
                for (const auto& T : res.poses) poses.emplace_back(T.matrix());
                for (const auto& st : res.steps) rmse.push_back(st.rmse), iterations.push_back(st.iterations);
            };
            if (scans.empty()) throw std::invalid_argument("no scans");
            if (scans[0].cols() == 2)
                run(std::integral_constant<int, 2>{});
            else
                run(std::integral_constant<int, 3>{});
            py::dict d;
            d["poses"] = poses;
            d["rmse"] = rmse;
            d["iterations"] = iterations;
            return d;
        },
        py::arg("scans"), py::arg("metric") = "point_to_plane", py::arg("max_iterations") = 50,
        py::arg("max_correspondence") = 1.0, py::arg("trim") = 1.0, py::arg("tolerance") = 1e-6,
        py::arg("normal_k") = 10, py::arg("voxel") = 0.0, py::arg("constant_velocity") = true,
        "scan-to-scan ICP chained (section 5.4); returns {poses (T_world_scan), rmse, iterations}");
}

}  // namespace sensor_fusion::python
