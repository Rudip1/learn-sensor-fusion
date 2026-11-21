#pragma once
/// @file registration.hpp
/// Rigid registration of point clouds: the closed-form alignment of matched points (SVD), ICP with the
/// point-to-point and point-to-plane metrics, and LiDAR odometry by chaining scan-to-scan ICP.
/// Templates over the dimension (2 for planar scans, 3 for 3-D LiDAR).
/// Theory: 1_theory/05_registration.md. Equation numbers refer to that file.

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "sensor_fusion/cloud.hpp"

namespace sensor_fusion::reg {

template <int Dim>
using Rigid = Eigen::Transform<double, Dim, Eigen::Isometry>;

/// Rigid transform T minimizing sum_i |dst_i - T src_i|^2 for matched columns (Arun / Umeyama, without
/// scale): SVD of the cross-covariance, with the reflection excluded. Eqs. (5.2)-(5.3).
template <int Dim>
Rigid<Dim> best_rigid(const cloud::Points<Dim>& src, const cloud::Points<Dim>& dst);

enum class Metric { PointToPoint, PointToPlane };

struct IcpOptions {
    Metric metric = Metric::PointToPoint;
    int max_iterations = 50;
    double max_correspondence = 1.0;  ///< pairs farther apart than this are rejected [m]
    double trim = 1.0;                ///< keep this fraction of the closest pairs (1 = all within the gate)
    double tolerance = 1e-6;          ///< stop when the update moves points by less than this [m, rad]
    int normal_k = 10;                ///< neighbours for the target normals (point-to-plane)
};

template <int Dim>
struct IcpResult {
    Rigid<Dim> T = Rigid<Dim>::Identity();  ///< maps source points into the target frame
    int iterations = 0;
    bool converged = false;
    double rmse = 0.0;             ///< RMS residual of the kept pairs at the end (point or plane distance)
    int pairs = 0;                 ///< pairs kept in the last iteration
    std::vector<double> history;   ///< rmse after every iteration
};

/// ICP (Algorithm 5.1): alternate nearest-neighbour matching in a k-d tree of the target and the
/// minimization of eq. (5.2) (point-to-point, closed form) or eq. (5.4) (point-to-plane, one Gauss-Newton
/// step of the linearized problem (5.5)), starting from `initial`.
template <int Dim>
IcpResult<Dim> icp(const cloud::Points<Dim>& source, const cloud::Points<Dim>& target, const Rigid<Dim>& initial,
                   const IcpOptions& options = {});

struct OdometryOptions {
    IcpOptions icp;
    double voxel = 0.0;              ///< downsample each scan with this leaf before matching (0 = no)
    bool constant_velocity = true;   ///< initial guess = previous relative motion (else identity)
};

template <int Dim>
struct OdometryResult {
    std::vector<Rigid<Dim>> poses;   ///< T_world_scan_k, first = identity (eq. 5.6)
    std::vector<IcpResult<Dim>> steps;
};

/// LiDAR odometry (section 5.4): ICP of every scan against the previous one, chained.
template <int Dim>
OdometryResult<Dim> lidar_odometry(const std::vector<cloud::Points<Dim>>& scans, const OdometryOptions& options = {});

}  // namespace sensor_fusion::reg
