// Chapter 8 bindings: ground segmentation, DBSCAN, persistence and the dynamic-object filter.
#include <pybind11/stl.h>

#include "bindings.hpp"
#include "sensor_fusion/objects.hpp"

namespace sensor_fusion::python {

namespace {
using RowsXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
Eigen::Matrix3Xd cols3(const RowsXd& r) {
    if (r.cols() != 3) throw std::invalid_argument("points must have shape (n, 3)");
    return r.transpose();
}
}  // namespace

void bind_objects(py::module_ m) {
    using namespace objects;
    m.def(
        "segment_ground",
        [](const RowsXd& pts, double threshold, int iterations, double max_tilt_deg, std::uint32_t seed) {
            const auto r = segment_ground(cols3(pts), threshold, iterations, max_tilt_deg * 0.017453292519943295, seed);
            py::dict d;
            d["normal"] = r.plane.n;
            d["d"] = r.plane.d;
            d["ground"] = r.ground;
            d["inliers"] = r.inliers;
            return d;
        },
        py::arg("points"), py::arg("threshold") = 0.15, py::arg("iterations") = 500, py::arg("max_tilt_deg") = 15.0,
        py::arg("seed") = 0, "RANSAC ground plane (Algorithm 8.1); returns {normal, d, ground, inliers}");
    m.def(
        "ground_by_height", [](const RowsXd& pts, double z_max) { return ground_by_height(cols3(pts), z_max); },
        py::arg("points"), py::arg("z_max"), "naive ground: z < z_max in the sensor frame");
    m.def(
        "dbscan",
        [](const RowsXd& pts, double eps, int min_points) {
            if (pts.cols() == 2) return dbscan<2>(cloud::Points<2>(pts.transpose()), eps, min_points);
            return dbscan<3>(cols3(pts), eps, min_points);
        },
        py::arg("points"), py::arg("eps"), py::arg("min_points"), "DBSCAN labels (Algorithm 8.2), -1 = noise");
    m.def(
        "persistence",
        [](const std::vector<RowsXd>& scans, double voxel, bool dilate) {
            std::vector<Eigen::Matrix3Xd> s;
            for (const auto& x : scans) s.push_back(cols3(x));
            return persistence(s, voxel, dilate);
        },
        py::arg("world_scans"), py::arg("voxel"), py::arg("dilate") = false, "per-point fraction of scans occupying its voxel, eq. (8.3)");
    m.def(
        "dynamic_points",
        [](const std::vector<RowsXd>& scans, const std::vector<std::vector<int>>& labels, double voxel, double ratio,
           bool dilate) {
            std::vector<Eigen::Matrix3Xd> s;
            for (const auto& x : scans) s.push_back(cols3(x));
            return dynamic_points(s, labels, voxel, ratio, dilate);
        },
        py::arg("world_scans"), py::arg("cluster_labels"), py::arg("voxel") = 0.5, py::arg("static_ratio") = 0.3,
        py::arg("dilate") = false,
        "per-point dynamic flags, eq. (8.4)");
}

}  // namespace sensor_fusion::python
