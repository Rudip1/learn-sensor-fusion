// Python module sensor_fusion._core: one submodule per topic.
#include <pybind11/pybind11.h>

#include "bindings.hpp"

PYBIND11_MODULE(_core, m) {
    m.doc() = "C++ core of the sensor_fusion learning module";
    sensor_fusion::python::bind_nmea(m.def_submodule("nmea", "NMEA 0183 decoding (chapter 1)"));
    sensor_fusion::python::bind_gnss(m.def_submodule("gnss", "Pseudorange positioning (chapter 1)"));
    sensor_fusion::python::bind_geo(m.def_submodule("geo", "WGS-84, ECEF, ENU and UTM (chapter 2)"));
    sensor_fusion::python::bind_cloud(m.def_submodule("cloud", "PCD files, k-d tree, voxel grid, normals (chapter 4)"));
    sensor_fusion::python::bind_sim(m.def_submodule("sim", "synthetic scene, ray-cast LiDAR and camera"));
}
