#pragma once
/// @file sim.hpp
/// A minimal synthetic world for experiments with exact ground truth: a flat ground plane and oriented
/// boxes (cars, buildings, posts), a ray-cast spinning LiDAR and a ray-cast pinhole camera.
/// The models are described in 1_theory/04_point_clouds.md (LiDAR) and 07_map_building.md (camera).

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cstdint>
#include <vector>

namespace sensor_fusion::sim {

/// An axis-aligned box rotated by `yaw` about the vertical axis.
struct Box {
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    Eigen::Vector3d half_extent = Eigen::Vector3d::Ones();
    double yaw = 0.0;
    Eigen::Vector3d color = Eigen::Vector3d(0.6, 0.6, 0.6);  ///< RGB in [0, 1]
    int id = 0;                                               ///< label written into simulated data
};

struct Scene {
    double ground_z = 0.0;
    Eigen::Vector3d ground_color = Eigen::Vector3d(0.35, 0.35, 0.38);
    double ground_checker = 2.0;  ///< edge of the ground checkerboard [m] (0 = plain); gives texture to images
    std::vector<Box> boxes;
};

/// A parking-lot-like scene: two rows of cars along a 60 m lane, buildings on both sides, light posts.
Scene demo_scene();

/// Spinning multi-beam LiDAR. A single beam at elevation 0 is a planar (2-D) scanner.
struct LidarModel {
    std::vector<double> elevations;  ///< beam elevations [rad]
    int azimuth_steps = 720;         ///< rays per revolution
    double min_range = 0.5;
    double max_range = 60.0;
    double range_noise = 0.0;        ///< 1-sigma Gaussian noise on the range [m]

    /// `beams` beams evenly spaced in [min_elev, max_elev] (radians).
    static LidarModel spinning(int beams, double min_elev, double max_elev, int azimuth_steps = 720,
                               double max_range = 60.0);
    /// One beam at elevation 0.
    static LidarModel planar(int azimuth_steps = 720, double max_range = 30.0);
};

struct Hit {
    double range = 0.0;
    int label = -2;  ///< box id, -1 for the ground, -2 for no hit
    Eigen::Vector3d normal = Eigen::Vector3d::Zero();  ///< surface normal at the hit (world frame)
};

/// First intersection of the ray o + t d (t > t_min) with the scene, eq. (4.1). d must be a unit vector.
Hit cast_ray(const Scene& scene, const Eigen::Vector3d& o, const Eigen::Vector3d& d, double t_min, double t_max);

struct LidarScan {
    Eigen::Matrix3Xd points;   ///< in the sensor frame
    Eigen::VectorXi labels;    ///< box id or -1 (ground)
    Eigen::Matrix3Xd normals;  ///< true surface normals, sensor frame
};

/// Simulate one revolution from the sensor pose T_world_sensor. Rays that hit nothing are dropped.
LidarScan simulate_lidar(const Scene& scene, const Eigen::Isometry3d& T_world_sensor, const LidarModel& model,
                         std::uint32_t seed = 0);

}  // namespace sensor_fusion::sim
