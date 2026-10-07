#include <pybind11/stl.h>

#include "bindings.hpp"
#include "sensor_fusion/sim.hpp"

namespace sensor_fusion::python {

Eigen::Isometry3d to_isometry(const Eigen::Matrix4d& T) {
    Eigen::Isometry3d iso = Eigen::Isometry3d::Identity();
    iso.linear() = T.topLeftCorner<3, 3>();
    iso.translation() = T.topRightCorner<3, 1>();
    return iso;
}

void bind_sim(py::module_ m) {
    using namespace sensor_fusion::sim;
    py::class_<Box>(m, "Box")
        .def(py::init([](const Eigen::Vector3d& center, const Eigen::Vector3d& half_extent, double yaw,
                         const Eigen::Vector3d& color, int id) {
                 return Box{center, half_extent, yaw, color, id};
             }),
             py::arg("center"), py::arg("half_extent"), py::arg("yaw") = 0.0,
             py::arg("color") = Eigen::Vector3d(0.6, 0.6, 0.6), py::arg("id") = 0)
        .def_readwrite("center", &Box::center)
        .def_readwrite("half_extent", &Box::half_extent)
        .def_readwrite("yaw", &Box::yaw)
        .def_readwrite("color", &Box::color)
        .def_readwrite("id", &Box::id)
        .def("__repr__", [](const Box& b) {
            return "<Box id " + std::to_string(b.id) + " at (" + std::to_string(b.center.x()) + ", " +
                   std::to_string(b.center.y()) + ", " + std::to_string(b.center.z()) + ")>";
        });
    py::class_<Scene>(m, "Scene", "ground plane + boxes; assign a whole list to .boxes to change it")
        .def(py::init<>())
        .def_readwrite("ground_z", &Scene::ground_z)
        .def_readwrite("ground_color", &Scene::ground_color)
        .def_readwrite("ground_checker", &Scene::ground_checker)
        .def_readwrite("boxes", &Scene::boxes);
    m.def("demo_scene", &demo_scene, "parking-lot-like scene: cars, buildings, posts");

    py::class_<LidarModel>(m, "LidarModel")
        .def(py::init<>())
        .def_readwrite("elevations", &LidarModel::elevations, "beam elevations [rad]")
        .def_readwrite("azimuth_steps", &LidarModel::azimuth_steps)
        .def_readwrite("min_range", &LidarModel::min_range)
        .def_readwrite("max_range", &LidarModel::max_range)
        .def_readwrite("range_noise", &LidarModel::range_noise)
        .def_static(
            "spinning",
            [](int beams, double min_deg, double max_deg, int steps, double max_range) {
                return LidarModel::spinning(beams, min_deg * 0.017453292519943295, max_deg * 0.017453292519943295,
                                            steps, max_range);
            },
            py::arg("beams"), py::arg("min_elevation_deg"), py::arg("max_elevation_deg"),
            py::arg("azimuth_steps") = 720, py::arg("max_range") = 60.0)
        .def_static("planar", &LidarModel::planar, py::arg("azimuth_steps") = 720, py::arg("max_range") = 30.0);

    py::class_<Hit>(m, "Hit")
        .def_readonly("range", &Hit::range)
        .def_readonly("label", &Hit::label)
        .def_readonly("normal", &Hit::normal);
    m.def("cast_ray", &cast_ray, py::arg("scene"), py::arg("origin"), py::arg("direction"), py::arg("t_min") = 0.0,
          py::arg("t_max") = 1e9, "first intersection of a ray with the scene, eq. (4.1)");
    m.def(
        "simulate_lidar",
        [](const Scene& scene, const Eigen::Matrix4d& T_world_sensor, const LidarModel& model, std::uint32_t seed) {
            const LidarScan s = simulate_lidar(scene, to_isometry(T_world_sensor), model, seed);
            py::dict out;
            out["points"] = RowsX3(s.points.transpose());
            out["labels"] = Eigen::VectorXi(s.labels);
            out["normals"] = RowsX3(s.normals.transpose());
            return out;
        },
        py::arg("scene"), py::arg("T_world_sensor"), py::arg("model"), py::arg("seed") = 0,
        "one revolution; returns {'points': (n,3) sensor frame, 'labels': (n,), 'normals': (n,3)}");
}

}  // namespace sensor_fusion::python
