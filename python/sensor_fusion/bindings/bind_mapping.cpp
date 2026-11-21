// Chapter 7 bindings: sensor poses, aggregation, projection and colouring. Images are (h, w, 3) float arrays
// as returned by sim.render_camera()["rgb"].
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "bindings.hpp"
#include "sensor_fusion/mapping.hpp"

namespace sensor_fusion::python {

namespace {
using RowsXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

Eigen::Matrix3Xd cols3(const RowsXd& r) {
    if (r.cols() != 3) throw std::invalid_argument("points must have shape (n, 3)");
    return r.transpose();
}

sim::RgbImage to_image(const py::array_t<float, py::array::c_style | py::array::forcecast>& a) {
    if (a.ndim() != 3 || a.shape(2) != 3) throw std::invalid_argument("image must have shape (h, w, 3)");
    sim::RgbImage img;
    img.height = int(a.shape(0));
    img.width = int(a.shape(1));
    img.rgb = Eigen::Map<const Eigen::Matrix<float, 3, Eigen::Dynamic>>(a.data(), 3, Eigen::Index(img.width) * img.height);
    return img;
}
}  // namespace

void bind_mapping(py::module_ m) {
    using namespace mapping;
    m.def(
        "sensor_pose",
        [](const Eigen::Vector3d& p, const Eigen::Matrix4d& T_vs) -> Eigen::Matrix4d {
            return sensor_pose(p, to_isometry(T_vs)).matrix();
        },
        py::arg("vehicle_xy_yaw"), py::arg("T_vehicle_sensor"), "T_world_sensor, eq. (7.1)");
    m.def(
        "aggregate",
        [](const std::vector<RowsXd>& scans, const std::vector<Eigen::Matrix4d>& poses, double voxel) {
            std::vector<Eigen::Matrix3Xd> s;
            std::vector<Eigen::Isometry3d> T;
            for (const auto& x : scans) s.push_back(cols3(x));
            for (const auto& x : poses) T.push_back(to_isometry(x));
            const Map r = aggregate(s, T, voxel);
            py::dict d;
            d["points"] = RowsX3(r.points.transpose());
            d["scan"] = Eigen::VectorXi(r.scan);
            return d;
        },
        py::arg("scans"), py::arg("poses"), py::arg("voxel") = 0.0, "stack scans in the world frame, eq. (7.2)");
    m.def(
        "cloud_to_cloud_rms", [](const RowsXd& a, const RowsXd& b) { return cloud_to_cloud_rms(cols3(a), cols3(b)); },
        py::arg("cloud"), py::arg("reference"));
    m.def(
        "project",
        [](const RowsXd& pts, const Eigen::Matrix4d& T_wc, const sim::Camera& cam) {
            const Projection p = project(cols3(pts), to_isometry(T_wc), cam);
            py::dict d;
            d["pixels"] = RowsX2(p.pixels.transpose());
            d["depth"] = Eigen::VectorXd(p.depth);
            d["inside"] = p.inside;
            return d;
        },
        py::arg("points"), py::arg("T_world_camera"), py::arg("camera"), "pinhole projection, eq. (7.3)");
    m.def(
        "colorize",
        [](const RowsXd& pts, const std::vector<py::array_t<float, py::array::c_style | py::array::forcecast>>& images,
           const std::vector<Eigen::Matrix4d>& poses, const sim::Camera& cam, bool occlusion_check, double tol, int window) {
            std::vector<sim::RgbImage> imgs;
            std::vector<Eigen::Isometry3d> T;
            for (const auto& a : images) imgs.push_back(to_image(a));
            for (const auto& x : poses) T.push_back(to_isometry(x));
            const Colors c = colorize(cols3(pts), imgs, T, cam, occlusion_check, tol, window);
            py::dict d;
            d["rgb"] = RowsX3(c.rgb.transpose());
            d["camera"] = c.camera;
            return d;
        },
        py::arg("points"), py::arg("images"), py::arg("poses"), py::arg("camera"), py::arg("occlusion_check") = true,
        py::arg("tolerance") = 0.3, py::arg("window") = 2, "Algorithm 7.1; returns {'rgb': (n, 3), 'camera': list}");
}

}  // namespace sensor_fusion::python
