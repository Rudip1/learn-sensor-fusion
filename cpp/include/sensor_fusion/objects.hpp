#pragma once
/// @file objects.hpp
/// Dynamic objects in LiDAR scans: ground removal (RANSAC plane, or a naive height threshold), Euclidean
/// clustering with DBSCAN, and a temporal-consistency filter that labels clusters as static or dynamic from
/// how persistently their voxels are occupied over a window of scans.
/// Theory: 1_theory/08_dynamic_objects.md. Equation numbers refer to that file.

#include <Eigen/Core>
#include <cstdint>
#include <vector>

#include "sensor_fusion/cloud.hpp"

namespace sensor_fusion::objects {

/// Plane n . p + d = 0 with unit normal n (n_z >= 0).
struct Plane {
    Eigen::Vector3d n = Eigen::Vector3d::UnitZ();
    double d = 0.0;
    double distance(const Eigen::Vector3d& p) const { return n.dot(p) + d; }
};

struct GroundResult {
    Plane plane;
    std::vector<bool> ground;  ///< per point
    int inliers = 0;
};

/// Ground plane by RANSAC (Algorithm 8.1): planes through 3 random points whose normal is within `max_tilt`
/// [rad] of the z axis, scored by the number of points within `threshold` (eq. 8.1); the best is refitted to its
/// inliers by PCA (chapter 4). Points within `threshold` of the final plane are ground.
GroundResult segment_ground(const Eigen::Matrix3Xd& points, double threshold = 0.15, int iterations = 500,
                            double max_tilt = 0.26, std::uint32_t seed = 0);

/// The naive alternative: ground = points lower than `z_max` in the sensor frame.
std::vector<bool> ground_by_height(const Eigen::Matrix3Xd& points, double z_max);

/// DBSCAN (Algorithm 8.2): points with at least `min_points` neighbours within `eps` (themselves included) are
/// core points; clusters are the connected components of core points plus the border points they reach.
/// Returns a label per point: 0..k-1, or -1 for noise. Eq. (8.2).
template <int Dim>
std::vector<int> dbscan(const cloud::Points<Dim>& points, double eps, int min_points);

struct Cluster {
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    Eigen::Vector3d min = Eigen::Vector3d::Zero(), max = Eigen::Vector3d::Zero();  ///< axis-aligned box
    std::vector<int> indices;
};

/// Clusters from DBSCAN labels (noise skipped), in label order.
std::vector<Cluster> make_clusters(const Eigen::Matrix3Xd& points, const std::vector<int>& labels);

/// Persistence of every point of every scan: the fraction of the scans in the window in which the voxel of the
/// point is occupied (with `dilate`, the most occupied of the voxel and its 26 neighbours). Scans must be in the
/// world frame. Eq. (8.3).
std::vector<Eigen::VectorXd> persistence(const std::vector<Eigen::Matrix3Xd>& world_scans, double voxel,
                                         bool dilate = false);

/// Temporal-consistency filter (eq. 8.4): a cluster is dynamic if the mean persistence of its points is below
/// `static_ratio`. Returns a per-point flag for every scan (noise and ground points are never dynamic).
std::vector<std::vector<bool>> dynamic_points(const std::vector<Eigen::Matrix3Xd>& world_scans,
                                              const std::vector<std::vector<int>>& cluster_labels, double voxel = 0.5,
                                              double static_ratio = 0.3, bool dilate = false);

}  // namespace sensor_fusion::objects
