#pragma once
/// @file mapping.hpp
/// Map building: aggregating scans into one cloud with the poses of a localization system, projecting the map
/// into camera images, and colouring it with an occlusion test.
/// Theory: 1_theory/07_map_building.md. Equation numbers refer to that file.

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "sensor_fusion/sim.hpp"

namespace sensor_fusion::mapping {

/// 3-D pose of a vehicle-mounted sensor from the planar vehicle pose (x, y, yaw) and the fixed mounting
/// T_vehicle_sensor: T_world_sensor = T_world_vehicle(x, y, yaw) T_vehicle_sensor. Eq. (7.1).
Eigen::Isometry3d sensor_pose(const Eigen::Vector3d& vehicle_xy_yaw, const Eigen::Isometry3d& T_vehicle_sensor);

struct Map {
    Eigen::Matrix3Xd points;  ///< world frame
    Eigen::VectorXi scan;     ///< index of the scan each point came from
};

/// Transform every scan into the world frame with its pose and stack them (eq. 7.2); if voxel > 0, downsample
/// the result (chapter 4), keeping for each voxel the scan of its first point.
Map aggregate(const std::vector<Eigen::Matrix3Xd>& scans, const std::vector<Eigen::Isometry3d>& T_world_sensor,
              double voxel = 0.0);

/// RMS distance from every point of `cloud` to its nearest neighbour in `reference` (a measure of map blur).
double cloud_to_cloud_rms(const Eigen::Matrix3Xd& cloud, const Eigen::Matrix3Xd& reference);

struct Projection {
    Eigen::Matrix2Xd pixels;    ///< (u, v), pixel coordinates (centre of the top-left pixel at (0.5, 0.5))
    Eigen::VectorXd depth;      ///< z in the camera frame
    std::vector<bool> inside;   ///< in front of the camera and inside the image
};

/// Pinhole projection of world points: X_c = T_world_camera^-1 X, (u, v) = (K X_c) / z. Eq. (7.3).
Projection project(const Eigen::Matrix3Xd& world, const Eigen::Isometry3d& T_world_camera, const sim::Camera& camera);

struct Colors {
    Eigen::Matrix3Xd rgb;     ///< colour of every point (0 where no camera sees it)
    std::vector<int> camera;  ///< index of the image that coloured the point, -1 if none
};

/// Colour map points from images (Algorithm 7.1). For each image, project all points; with `occlusion_check`,
/// a point counts as visible only if its depth is within `tolerance` of the smallest depth projected into a
/// window of (2 r + 1)^2 pixels around it (a z-buffer of the map itself). Each point takes the colour from the
/// image in which it is closest to the camera.
Colors colorize(const Eigen::Matrix3Xd& world, const std::vector<sim::RgbImage>& images,
                const std::vector<Eigen::Isometry3d>& T_world_camera, const sim::Camera& camera,
                bool occlusion_check = true, double tolerance = 0.3, int window = 2);

}  // namespace sensor_fusion::mapping
