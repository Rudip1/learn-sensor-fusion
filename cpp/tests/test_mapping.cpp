// Chapter 7 tests: sensor poses from vehicle poses, aggregation, the displacement a heading error causes, pinhole
// projection against the ray-cast camera, and colouring with and without the occlusion test.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>

#include "sensor_fusion/mapping.hpp"
#include "sensor_fusion/sim.hpp"

using namespace sensor_fusion;
using Catch::Approx;

namespace {

Eigen::Isometry3d mount(double x, double y, double z) {
    Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
    T.translation() = Eigen::Vector3d(x, y, z);
    return T;
}

/// Camera looking along the vehicle's x axis (camera z forward, x right = -vehicle y, y down = -vehicle z).
Eigen::Isometry3d forward_camera(double height) {
    Eigen::Isometry3d T = mount(0, 0, height);
    T.linear() << 0, 0, 1, -1, 0, 0, 0, -1, 0;
    return T;
}

sim::Camera camera() {
    sim::Camera c;
    c.K << 300, 0, 200, 0, 300, 150, 0, 0, 1;
    c.width = 400;
    c.height = 300;
    return c;
}

}  // namespace

TEST_CASE("sensor pose and aggregation", "[mapping][ch7]") {
    const Eigen::Isometry3d T = mapping::sensor_pose(Eigen::Vector3d(10, 5, M_PI / 2), mount(1.0, 0, 1.8));
    CHECK((T.translation() - Eigen::Vector3d(10, 6, 1.8)).norm() < 1e-12);  // the mount is rotated with the vehicle
    Eigen::Matrix3Xd a(3, 2), b(3, 1);
    a << 1, 2, 0, 0, 0, 0;
    b << 0, 1, 0;
    const auto m = mapping::aggregate({a, b}, {mount(0, 0, 0), T});
    REQUIRE(m.points.cols() == 3);
    CHECK(m.scan(2) == 1);
    CHECK((m.points.col(2) - T * Eigen::Vector3d(0, 1, 0)).norm() < 1e-12);
}

TEST_CASE("a heading error displaces a point at range r by 2 r sin(dpsi / 2)", "[mapping][ch7]") {
    Eigen::Matrix3Xd p(3, 3);
    p << 5, 20, 40, 0, 0, 0, 0, 0, 0;
    const double d = 1.0 * M_PI / 180;
    const auto good = mapping::aggregate({p}, {mapping::sensor_pose(Eigen::Vector3d(0, 0, 0), mount(0, 0, 0))});
    const auto bad = mapping::aggregate({p}, {mapping::sensor_pose(Eigen::Vector3d(0, 0, d), mount(0, 0, 0))});
    for (int i = 0; i < 3; ++i)
        CHECK((bad.points.col(i) - good.points.col(i)).norm() == Approx(2 * p(0, i) * std::sin(d / 2)));
    CHECK((bad.points.col(2) - good.points.col(2)).norm() == Approx(0.698).epsilon(1e-3));  // 1 degree at 40 m
}

TEST_CASE("projection inverts the ray-cast camera", "[mapping][ch7]") {
    const sim::Scene scene = sim::demo_scene();
    const sim::Camera cam = camera();
    const Eigen::Isometry3d Twc = mapping::sensor_pose(Eigen::Vector3d(-20, 0, 0.1), forward_camera(1.5));
    const sim::RgbImage img = sim::render_camera(scene, Twc, cam);
    // back-project a few pixel centres with the rendered depth and project them again
    Eigen::Matrix3Xd X(3, 4);
    const int pix[4][2] = {{200, 160}, {20, 280}, {390, 220}, {123, 201}};  // below the horizon (v = 150)
    for (int k = 0; k < 4; ++k) {
        const int u = pix[k][0], v = pix[k][1];
        const double z = img.depth(Eigen::Index(v) * cam.width + u);
        REQUIRE(std::isfinite(z));
        X.col(k) = Twc * (z * cam.K.inverse() * Eigen::Vector3d(u + 0.5, v + 0.5, 1.0));
    }
    const auto p = mapping::project(X, Twc, cam);
    for (int k = 0; k < 4; ++k) {
        CHECK(p.inside[std::size_t(k)]);
        CHECK((p.pixels.col(k) - Eigen::Vector2d(pix[k][0] + 0.5, pix[k][1] + 0.5)).norm() < 1e-6);
    }
    // a point behind the camera is not inside
    const auto q = mapping::project(Eigen::Matrix3Xd(Twc * Eigen::Vector3d(0, 0, -5)), Twc, cam);
    CHECK_FALSE(q.inside[0]);
}

TEST_CASE("colouring needs the occlusion test", "[mapping][ch7]") {
    // a red box 6 m in front of a large blue wall; the LiDAR (offset sideways) sees parts of the wall that the
    // camera cannot see behind the box
    sim::Scene scene;
    scene.ground_checker = 0;
    sim::Box box, wall;
    box.center = {6, 0, 1};
    box.half_extent = {0.5, 1.5, 1};
    box.color = {0.9, 0.1, 0.1};
    box.id = 1;
    wall.center = {14, 0, 3};
    wall.half_extent = {0.5, 10, 3};
    wall.color = {0.1, 0.2, 0.9};
    wall.id = 2;
    scene.boxes = {box, wall};
    auto model = sim::LidarModel::spinning(32, -20 * M_PI / 180, 20 * M_PI / 180, 1440, 40.0);
    const Eigen::Isometry3d Tl = mount(0, 3, 1.5);
    const auto scan = sim::simulate_lidar(scene, Tl, model);
    const auto m = mapping::aggregate({scan.points}, {Tl});
    const sim::Camera cam = camera();
    const Eigen::Isometry3d Tc = forward_camera(1.5);
    const sim::RgbImage img = sim::render_camera(scene, Tc, cam);
    auto accuracy = [&](bool check) {
        const auto c = mapping::colorize(m.points, {img}, {Tc}, cam, check);
        const auto p = mapping::project(m.points, Tc, cam);
        int coloured = 0, correct = 0;
        for (Eigen::Index i = 0; i < m.points.cols(); ++i) {
            if (c.camera[std::size_t(i)] < 0) continue;
            ++coloured;
            const int u = int(p.pixels(0, i)), v = int(p.pixels(1, i));
            correct += img.label(Eigen::Index(v) * cam.width + u) == scan.labels(i);
        }
        return double(correct) / coloured;
    };
    CHECK(accuracy(true) > 0.97);
    CHECK(accuracy(false) < 0.9);
}
