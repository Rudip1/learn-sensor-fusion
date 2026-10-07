#include "sensor_fusion/sim.hpp"

#include <cmath>
#include <limits>
#include <random>

namespace sensor_fusion::sim {

namespace {
constexpr double kPi = 3.14159265358979323846;
}

Scene demo_scene() {
    Scene s;
    int id = 1;
    auto add = [&](double x, double y, double z, double hx, double hy, double hz, double yaw, Eigen::Vector3d c) {
        Box b;
        b.center = {x, y, z};
        b.half_extent = {hx, hy, hz};
        b.yaw = yaw;
        b.color = c;
        b.id = id++;
        s.boxes.push_back(b);
    };
    const Eigen::Vector3d palette[] = {{0.80, 0.15, 0.15}, {0.15, 0.35, 0.80}, {0.90, 0.85, 0.80},
                                       {0.10, 0.10, 0.10}, {0.20, 0.60, 0.30}, {0.85, 0.60, 0.10}};
    // two rows of parked cars (4.5 x 1.8 x 1.5 m), parked across the lane
    for (int i = 0; i < 10; ++i) {
        const double x = -27.0 + 6.0 * i;
        if (i != 3) add(x, 6.0, 0.75, 0.9, 2.25, 0.75, 0.05 * std::sin(i), palette[i % 6]);
        if (i != 6) add(x + 1.0, -6.0, 0.75, 0.9, 2.25, 0.75, -0.04 * std::cos(i), palette[(i + 3) % 6]);
    }
    // buildings
    add(0.0, 18.0, 5.0, 30.0, 3.0, 5.0, 0.0, {0.75, 0.70, 0.62});
    add(-10.0, -19.0, 4.0, 18.0, 3.0, 4.0, 0.0, {0.62, 0.66, 0.72});
    add(22.0, -17.0, 7.0, 6.0, 4.0, 7.0, 0.3, {0.70, 0.55, 0.50});
    // light posts
    for (int i = 0; i < 4; ++i) add(-24.0 + 16.0 * i, 10.5, 2.5, 0.12, 0.12, 2.5, 0.0, {0.3, 0.3, 0.3});
    return s;
}

LidarModel LidarModel::spinning(int beams, double min_elev, double max_elev, int azimuth_steps, double max_range) {
    LidarModel m;
    m.azimuth_steps = azimuth_steps;
    m.max_range = max_range;
    for (int i = 0; i < beams; ++i)
        m.elevations.push_back(beams == 1 ? 0.5 * (min_elev + max_elev)
                                          : min_elev + (max_elev - min_elev) * i / (beams - 1));
    return m;
}

LidarModel LidarModel::planar(int azimuth_steps, double max_range) {
    LidarModel m;
    m.azimuth_steps = azimuth_steps;
    m.max_range = max_range;
    m.elevations = {0.0};
    return m;
}

Hit cast_ray(const Scene& scene, const Eigen::Vector3d& o, const Eigen::Vector3d& d, double t_min, double t_max) {
    Hit best;
    double t_best = t_max;
    // ground plane z = ground_z
    if (d.z() < -1e-12) {
        const double t = (scene.ground_z - o.z()) / d.z();
        if (t > t_min && t < t_best) {
            t_best = t;
            best.label = -1;
            best.normal = Eigen::Vector3d::UnitZ();
        }
    }
    // boxes: slab test in the box frame (eq. 4.1)
    for (const Box& b : scene.boxes) {
        const double c = std::cos(b.yaw), s = std::sin(b.yaw);
        const Eigen::Vector3d po = o - b.center;
        const Eigen::Vector3d ol(c * po.x() + s * po.y(), -s * po.x() + c * po.y(), po.z());
        const Eigen::Vector3d dl(c * d.x() + s * d.y(), -s * d.x() + c * d.y(), d.z());
        double t0 = -std::numeric_limits<double>::infinity(), t1 = std::numeric_limits<double>::infinity();
        int axis0 = -1;
        double sign0 = 0.0;
        bool miss = false;
        for (int a = 0; a < 3; ++a) {
            if (std::abs(dl(a)) < 1e-12) {
                if (std::abs(ol(a)) > b.half_extent(a)) miss = true;
                continue;
            }
            double ta = (-b.half_extent(a) - ol(a)) / dl(a);
            double tb = (b.half_extent(a) - ol(a)) / dl(a);
            double sgn = -1.0;  // entering through the face at -half_extent
            if (ta > tb) {
                std::swap(ta, tb);
                sgn = 1.0;
            }
            if (ta > t0) {
                t0 = ta;
                axis0 = a;
                sign0 = sgn;
            }
            t1 = std::min(t1, tb);
        }
        if (miss || t0 > t1 || axis0 < 0) continue;
        if (t0 > t_min && t0 < t_best) {
            t_best = t0;
            best.label = b.id;
            Eigen::Vector3d nl = Eigen::Vector3d::Zero();
            nl(axis0) = sign0;
            best.normal = Eigen::Vector3d(c * nl.x() - s * nl.y(), s * nl.x() + c * nl.y(), nl.z());
        }
    }
    if (best.label != -2) best.range = t_best;
    return best;
}

LidarScan simulate_lidar(const Scene& scene, const Eigen::Isometry3d& T_world_sensor, const LidarModel& model,
                         std::uint32_t seed) {
    std::mt19937 rng(seed);
    std::normal_distribution<double> noise(0.0, model.range_noise > 0 ? model.range_noise : 1.0);
    const Eigen::Matrix3d R = T_world_sensor.linear();
    const Eigen::Vector3d o = T_world_sensor.translation();

    std::vector<Eigen::Vector3d> pts, nrm;
    std::vector<int> labels;
    for (double el : model.elevations) {
        for (int k = 0; k < model.azimuth_steps; ++k) {
            const double az = 2.0 * kPi * k / model.azimuth_steps;
            const Eigen::Vector3d ds(std::cos(el) * std::cos(az), std::cos(el) * std::sin(az), std::sin(el));
            const Hit h = cast_ray(scene, o, R * ds, model.min_range, model.max_range);
            if (h.label == -2) continue;
            const double r = h.range + (model.range_noise > 0 ? noise(rng) : 0.0);
            pts.push_back(r * ds);
            nrm.push_back(R.transpose() * h.normal);
            labels.push_back(h.label);
        }
    }
    LidarScan scan;
    scan.points.resize(3, static_cast<Eigen::Index>(pts.size()));
    scan.normals.resize(3, static_cast<Eigen::Index>(pts.size()));
    scan.labels.resize(static_cast<Eigen::Index>(pts.size()));
    for (std::size_t i = 0; i < pts.size(); ++i) {
        scan.points.col(static_cast<Eigen::Index>(i)) = pts[i];
        scan.normals.col(static_cast<Eigen::Index>(i)) = nrm[i];
        scan.labels(static_cast<Eigen::Index>(i)) = labels[i];
    }
    return scan;
}

}  // namespace sensor_fusion::sim
