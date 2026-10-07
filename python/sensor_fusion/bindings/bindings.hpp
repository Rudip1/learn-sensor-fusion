#pragma once
// Declarations shared by the binding translation units. Convention: the C++ library stores point sets as
// 3 x N (columns); Python sees N x 3 (rows), the layout numpy, scipy and Open3D use.
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace sensor_fusion::python {

namespace py = pybind11;

using RowsX3 = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using RowsX2 = Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor>;

void bind_nmea(py::module_ m);
void bind_gnss(py::module_ m);
void bind_geo(py::module_ m);
void bind_cloud(py::module_ m);
void bind_sim(py::module_ m);

/// 4x4 homogeneous matrix -> isometry (rotation part assumed orthonormal).
Eigen::Isometry3d to_isometry(const Eigen::Matrix4d& T);

}  // namespace sensor_fusion::python
