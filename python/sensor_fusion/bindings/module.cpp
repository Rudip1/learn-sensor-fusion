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
    sensor_fusion::python::bind_reg(m.def_submodule("reg", "registration: ICP and LiDAR odometry (chapter 5)"));
    sensor_fusion::python::bind_fusion(m.def_submodule("fusion", "GNSS + odometry EKF (chapter 6)"));
    sensor_fusion::python::bind_mapping(m.def_submodule("mapping", "aggregation, projection, colouring (chapter 7)"));
    sensor_fusion::python::bind_objects(m.def_submodule("objects", "ground removal, DBSCAN, dynamic objects (chapter 8)"));
}
