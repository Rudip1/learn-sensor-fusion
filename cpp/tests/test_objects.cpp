// Chapter 8 tests: RANSAC ground plane against a known tilted plane and against a height threshold on a pitched
// sensor, DBSCAN against a brute-force reference and on chains of points, persistence, and the whole
// dynamic-object pipeline on a simulated scene with a moving car and a pedestrian.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cmath>
#include <map>
#include <random>

#include "sensor_fusion/mapping.hpp"
#include "sensor_fusion/objects.hpp"
#include "sensor_fusion/sim.hpp"

using namespace sensor_fusion;
using Catch::Approx;

TEST_CASE("RANSAC finds a tilted ground plane among obstacles and ignores walls", "[objects][ch8]") {
    std::mt19937 gen(1);
    std::uniform_real_distribution<double> U(-20, 20);
    std::normal_distribution<double> N(0, 0.02);
    const Eigen::Vector3d n = Eigen::Vector3d(0.05, -0.03, 1).normalized();
    std::vector<Eigen::Vector3d> pts;
    for (int i = 0; i < 2000; ++i) {  // ground: n . p = -1.8
        const double x = U(gen), y = U(gen);
        pts.emplace_back(x, y, (-1.8 - n.x() * x - n.y() * y) / n.z() + N(gen));
    }
    for (int i = 0; i < 2500; ++i) pts.emplace_back(8 + N(gen), U(gen), U(gen) / 4);  // a large wall x = 8
    Eigen::Matrix3Xd P(3, Eigen::Index(pts.size()));
    for (std::size_t i = 0; i < pts.size(); ++i) P.col(Eigen::Index(i)) = pts[i];
    const auto r = objects::segment_ground(P);
    CHECK(std::acos(std::min(1.0, r.plane.n.dot(n))) < 0.2 * M_PI / 180);
    CHECK(r.plane.d == Approx(1.8).margin(0.01));
    CHECK(r.inliers >= 1990);
    CHECK(r.inliers < 2100);  // the wall (more points than the ground!) is not taken as ground
}

TEST_CASE("a height threshold fails on a pitched sensor; RANSAC does not", "[objects][ch8]") {
    const sim::Scene scene = sim::demo_scene();
    Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
    T.linear() = Eigen::AngleAxisd(3 * M_PI / 180, Eigen::Vector3d::UnitY()).toRotationMatrix();  // nose down 3 deg
    T.translation() = Eigen::Vector3d(-10, 0, 1.8);
    auto model = sim::LidarModel::spinning(16, -15 * M_PI / 180, 15 * M_PI / 180, 720, 60.0);
    model.range_noise = 0.02;
    const auto s = sim::simulate_lidar(scene, T, model);
    const auto ransac = objects::segment_ground(s.points);
    const auto height = objects::ground_by_height(s.points, -1.8 + 0.15);
    int ok_r = 0, ok_h = 0;
    for (Eigen::Index i = 0; i < s.points.cols(); ++i) {
        const bool truth = s.labels(i) == -1;
        ok_r += ransac.ground[std::size_t(i)] == truth;
        ok_h += height[std::size_t(i)] == truth;
    }
    const double n = double(s.points.cols());
    CHECK(ok_r / n > 0.95);
    CHECK(ok_h / n < 0.9);
}

namespace {
/// O(n^2) DBSCAN straight from the definition, the reference for the k-d tree version.
std::vector<int> brute_dbscan(const cloud::Points<2>& P, double eps, int min_points) {
    const int n = int(P.cols());
    std::vector<std::vector<int>> nb(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if ((P.col(i) - P.col(j)).norm() <= eps) nb[std::size_t(i)].push_back(j);
    std::vector<int> label(std::size_t(n), -1);
    int c = 0;
    for (int i = 0; i < n; ++i) {
        if (label[std::size_t(i)] >= 0 || int(nb[std::size_t(i)].size()) < min_points) continue;
        std::vector<int> stack{i};
        label[std::size_t(i)] = c;
        while (!stack.empty()) {
            const int j = stack.back();
            stack.pop_back();
            if (int(nb[std::size_t(j)].size()) < min_points) continue;  // border point: do not expand
            for (int q : nb[std::size_t(j)])
                if (label[std::size_t(q)] < 0) label[std::size_t(q)] = c, stack.push_back(q);
        }
        ++c;
    }
    return label;
}
}  // namespace

TEST_CASE("DBSCAN equals the definition; chains connect only below eps", "[objects][ch8]") {
    std::mt19937 gen(2);
    std::normal_distribution<double> N(0, 0.3);
    std::uniform_real_distribution<double> U(-10, 10);
    cloud::Points<2> P(2, 330);
    for (int i = 0; i < 150; ++i) P.col(i) << -4 + N(gen), N(gen);
    for (int i = 150; i < 300; ++i) P.col(i) << 4 + N(gen), 1 + N(gen);
    for (int i = 300; i < 330; ++i) P.col(i) << U(gen), 6 + U(gen);  // sparse noise
    const auto a = objects::dbscan<2>(P, 0.5, 4);
    const auto b = brute_dbscan(P, 0.5, 4);
    // same partition up to the numbering of the clusters
    std::map<int, int> map;
    for (std::size_t i = 0; i < a.size(); ++i) {
        CHECK((a[i] < 0) == (b[i] < 0));
        if (a[i] >= 0) {
            const auto [it, inserted] = map.emplace(a[i], b[i]);
            CHECK(it->second == b[i]);
        }
    }
    CHECK(*std::max_element(a.begin(), a.end()) == 1);  // two clusters
    // a chain of points 0.4 m apart is one cluster with eps 0.5; 0.6 m apart it is all noise
    for (double step : {0.4, 0.6}) {
        cloud::Points<2> C(2, 20);
        for (int i = 0; i < 20; ++i) C.col(i) << step * i, 0;
        const auto l = objects::dbscan<2>(C, 0.5, 3);
        if (step < 0.5)
            CHECK(std::all_of(l.begin(), l.end(), [](int v) { return v == 0; }));
        else
            CHECK(std::all_of(l.begin(), l.end(), [](int v) { return v == -1; }));
    }
}

TEST_CASE("persistence: 1 for what every scan sees, 1/N for what one scan sees", "[objects][ch8]") {
    std::vector<Eigen::Matrix3Xd> scans;
    for (int k = 0; k < 10; ++k) {
        Eigen::Matrix3Xd s(3, 2);
        s.col(0) << 1.1, 2.2, 0.3;            // a static point
        s.col(1) << 5.0 + k, 0.2, 0.3;        // a point moving 1 m per scan
        scans.push_back(s);
    }
    const auto p = objects::persistence(scans, 0.5);
    for (const auto& v : p) {
        CHECK(v(0) == Approx(1.0));
        CHECK(v(1) == Approx(0.1));
    }
}

TEST_CASE("dynamic objects in a simulated scene: a moving car and a pedestrian", "[objects][ch8]") {
    auto model = sim::LidarModel::spinning(16, -15 * M_PI / 180, 15 * M_PI / 180, 720, 60.0);
    model.range_noise = 0.02;
    Eigen::Isometry3d mount = Eigen::Isometry3d::Identity();
    mount.translation() = Eigen::Vector3d(0, 0, 1.8);
    std::vector<Eigen::Matrix3Xd> world;
    std::vector<std::vector<int>> labels;
    std::vector<Eigen::VectorXi> truth;
    for (int k = 0; k < 30; ++k) {  // 3 s at 10 Hz, ego vehicle standing at x = -10
        const double t = 0.1 * k;
        sim::Scene scene = sim::demo_scene();
        sim::Box car, ped;
        car.center = {15 - 5 * t, 2.5, 0.75}, car.half_extent = {2.25, 0.9, 0.75}, car.id = 100;
        ped.center = {8, -4 + 1.2 * t, 0.9}, ped.half_extent = {0.3, 0.3, 0.9}, ped.id = 101;
        scene.boxes.push_back(car);
        scene.boxes.push_back(ped);
        const Eigen::Isometry3d T = mapping::sensor_pose(Eigen::Vector3d(-10, 0, 0), mount);
        const auto s = sim::simulate_lidar(scene, T, model, std::uint32_t(k));
        const auto g = objects::segment_ground(s.points, 0.15, 500, 0.26, std::uint32_t(k));
        std::vector<int> idx;
        for (Eigen::Index i = 0; i < s.points.cols(); ++i)
            if (!g.ground[std::size_t(i)]) idx.push_back(int(i));
        cloud::Points<3> ng(3, Eigen::Index(idx.size()));
        for (std::size_t i = 0; i < idx.size(); ++i) ng.col(Eigen::Index(i)) = s.points.col(idx[i]);
        const auto l = objects::dbscan<3>(ng, 0.5, 3);
        std::vector<int> lab(std::size_t(s.points.cols()), -1);
        for (std::size_t i = 0; i < idx.size(); ++i) lab[std::size_t(idx[i])] = l[i];
        world.push_back(T * s.points);
        labels.push_back(lab);
        truth.push_back(s.labels);
    }
    const auto dyn = objects::dynamic_points(world, labels, 0.5, 0.3);
    int tp = 0, flagged = 0, moving = 0;
    for (std::size_t k = 0; k < dyn.size(); ++k)
        for (std::size_t i = 0; i < dyn[k].size(); ++i) {
            const bool m = truth[k](Eigen::Index(i)) >= 100;
            moving += m;
            flagged += dyn[k][i];
            tp += m && dyn[k][i];
        }
    CHECK(double(tp) / flagged > 0.95);  // precision
    CHECK(double(tp) / moving > 0.65);   // recall: part of the car's long sides looks static (section 8.4)
}
