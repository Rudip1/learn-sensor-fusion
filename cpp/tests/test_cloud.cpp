// Point-cloud basics (1_theory/04_point_clouds.md): k-d tree against brute force, voxel grid against a
// hand-worked example, PCA normals against known planes and the simulator's true normals, PCD round trips.
#include <Eigen/Dense>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>

#include "sensor_fusion/cloud.hpp"
#include "sensor_fusion/sim.hpp"

using Catch::Approx;
namespace cloud = sensor_fusion::cloud;
namespace sim = sensor_fusion::sim;

namespace {
template <int Dim>
cloud::Points<Dim> random_points(int n, unsigned seed, double scale = 10.0) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-scale, scale);
    cloud::Points<Dim> p(Dim, n);
    for (int i = 0; i < n; ++i)
        for (int d = 0; d < Dim; ++d) p(d, i) = u(rng);
    return p;
}
}  // namespace

TEST_CASE("k-d tree k-NN equals brute force (2-D and 3-D)", "[kdtree][ch4]") {
    for (int leaf : {1, 4, 16}) {
        const auto p3 = random_points<3>(2000, 1);
        const cloud::KdTree<3> t3(p3, leaf);
        const auto p2 = random_points<2>(2000, 2);
        const cloud::KdTree<2> t2(p2, leaf);
        const auto q3 = random_points<3>(100, 3, 12.0);
        const auto q2 = random_points<2>(100, 4, 12.0);
        std::vector<int> ia, ib;
        std::vector<double> da, db;
        for (int i = 0; i < 100; ++i) {
            for (int k : {1, 7, 30}) {
                t3.knn(q3.col(i), k, ia, da);
                cloud::brute_force_knn<3>(p3, q3.col(i), k, ib, db);
                REQUIRE(da.size() == db.size());
                for (std::size_t j = 0; j < da.size(); ++j) CHECK(da[j] == Approx(db[j]));
                t2.knn(q2.col(i), k, ia, da);
                cloud::brute_force_knn<2>(p2, q2.col(i), k, ib, db);
                for (std::size_t j = 0; j < da.size(); ++j) CHECK(da[j] == Approx(db[j]));
            }
        }
    }
}

TEST_CASE("k-d tree radius search equals brute force", "[kdtree][ch4]") {
    const auto p = random_points<3>(3000, 5);
    const cloud::KdTree<3> tree(p);
    const auto q = random_points<3>(50, 6);
    std::vector<int> idx;
    std::vector<double> d2;
    for (int i = 0; i < 50; ++i) {
        tree.radius(q.col(i), 2.5, idx, d2);
        int expected = 0;
        for (int j = 0; j < p.cols(); ++j) expected += (p.col(j) - q.col(i)).norm() <= 2.5;
        CHECK(static_cast<int>(idx.size()) == expected);
        CHECK(std::is_sorted(d2.begin(), d2.end()));
    }
}

TEST_CASE("k-d tree edge cases", "[kdtree][ch4]") {
    const cloud::KdTree<2> empty(cloud::Points<2>(2, 0));
    CHECK(empty.nearest(cloud::Point<2>(0, 0)).first == -1);
    // many identical points: must not recurse forever and must return them all
    cloud::Points<2> same = cloud::Points<2>::Ones(2, 100);
    const cloud::KdTree<2> t(same, 4);
    std::vector<int> idx;
    std::vector<double> d2;
    t.knn(cloud::Point<2>(1, 1), 50, idx, d2);
    CHECK(idx.size() == 50);
    CHECK(d2.back() == 0.0);
    t.knn(cloud::Point<2>(1, 1), 500, idx, d2);  // k larger than the tree
    CHECK(idx.size() == 100);
    // a query point that is in the cloud is its own nearest neighbour
    const auto p = random_points<3>(500, 9);
    const cloud::KdTree<3> t3(p);
    CHECK(t3.nearest(p.col(123)).first == 123);
    CHECK(t3.depth() <= 10);  // balanced: ceil(log2(500 / 10)) + 1 levels at most
}

TEST_CASE("voxel grid: hand-worked example of eqs. (4.2)-(4.3)", "[voxel][ch4]") {
    cloud::Points<2> p(2, 5);
    // leaf 1: (0.2,0.2) and (0.8,0.6) share voxel (0,0); (1.5,0.5) is in (1,0); (-0.5,-0.5) in (-1,-1);
    // (-0.1, 0.9) in (-1, 0)
    p << 0.2, 0.8, 1.5, -0.5, -0.1,
         0.2, 0.6, 0.5, -0.5, 0.9;
    const auto r = cloud::voxel_downsample<2>(p, 1.0);
    REQUIRE(r.centroids.cols() == 4);
    CHECK(r.centroids.col(0).isApprox(cloud::Point<2>(0.5, 0.4)));
    CHECK(r.count[0] == 2);
    CHECK(r.keys[2] == Eigen::Matrix<long long, 2, 1>(-1, -1));
    CHECK(r.keys[3] == Eigen::Matrix<long long, 2, 1>(-1, 0));
    CHECK(r.voxel_of_point == std::vector<int>{0, 0, 1, 2, 3});
    // shifting the origin by half a voxel regroups the points
    const auto s = cloud::voxel_downsample<2>(p, 1.0, cloud::Point<2>(0.5, 0.5));
    CHECK(s.keys[0] == Eigen::Matrix<long long, 2, 1>(-1, -1));
    CHECK(s.voxel_of_point[0] != s.voxel_of_point[1]);
}

TEST_CASE("voxel grid keeps one point per voxel and preserves the centroid of the cloud", "[voxel][ch4]") {
    const auto p = random_points<3>(5000, 10, 5.0);
    const auto r = cloud::voxel_downsample<3>(p, 0.7);
    Eigen::Vector3d weighted = Eigen::Vector3d::Zero();
    for (int v = 0; v < r.centroids.cols(); ++v) weighted += r.count[static_cast<std::size_t>(v)] * r.centroids.col(v);
    CHECK((weighted / p.cols() - p.rowwise().mean()).norm() < 1e-12);
    // no two centroids share a voxel key
    for (std::size_t i = 0; i < r.keys.size(); ++i)
        for (std::size_t j = i + 1; j < r.keys.size(); ++j) CHECK(r.keys[i] != r.keys[j]);
}

TEST_CASE("PCA normal of noiseless planar points is the plane normal", "[normals][ch4]") {
    std::mt19937 rng(4);
    std::uniform_real_distribution<double> u(-1, 1);
    const Eigen::Vector3d n = Eigen::Vector3d(0.3, -0.5, 0.8).normalized();
    const Eigen::Vector3d a = n.unitOrthogonal(), b = n.cross(a);
    cloud::Points<3> p(3, 200);
    for (int i = 0; i < 200; ++i) p.col(i) = Eigen::Vector3d(1, 2, 3) + 4 * u(rng) * a + 2 * u(rng) * b;
    const auto f = cloud::fit_plane<3>(p);
    CHECK(std::abs(f.normal.dot(n)) == Approx(1.0).epsilon(1e-12));
    CHECK(f.variation < 1e-20);
    // a 2-D line
    cloud::Points<2> l(2, 50);
    for (int i = 0; i < 50; ++i) l.col(i) = cloud::Point<2>(i * 0.1, 2.0 + 0.5 * i * 0.1);
    const auto fl = cloud::fit_plane<2>(l);
    CHECK(std::abs(fl.normal.dot(cloud::Point<2>(-0.5, 1.0).normalized())) == Approx(1.0));
}

TEST_CASE("estimated normals agree with the simulator's true normals and face the sensor", "[normals][ch4]") {
    const sim::Scene scene = sim::demo_scene();
    Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
    T.translation() = Eigen::Vector3d(0, 0, 1.8);
    const auto scan = sim::simulate_lidar(scene, T, sim::LidarModel::spinning(32, -0.4, 0.2, 1024));
    const cloud::KdTree<3> tree(scan.points);
    const auto nr = cloud::estimate_normals<3>(scan.points, tree, 12);
    int flat = 0, good = 0;
    for (int i = 0; i < scan.points.cols(); ++i) {
        CHECK(nr.normals.col(i).dot(-scan.points.col(i)) >= 0.0);  // oriented towards the sensor at the origin
        if (nr.variation(i) < 1e-6) {  // neighbourhood lies on one face
            // ...unless all k neighbours lie on one line (a single scan ring far away): then two eigenvalues
            // vanish and the normal is undetermined within a plane (section 4.5, degenerate neighbourhoods)
            std::vector<int> idx;
            std::vector<double> d2;
            tree.knn(scan.points.col(i), 12, idx, d2);
            cloud::Points<3> nb(3, 12);
            for (int j = 0; j < 12; ++j) nb.col(j) = scan.points.col(idx[static_cast<std::size_t>(j)]);
            const Eigen::Vector3d ev =
                Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(cloud::mean_and_covariance<3>(nb).second).eigenvalues();
            if (ev(1) < 1e-6 * ev(2)) continue;
            // ...or mixes two surfaces (11 ground points and one at the foot of a car still look "flat")
            bool one_surface = true;
            for (int j : idx) one_surface &= scan.labels(j) == scan.labels(i);
            if (!one_surface) continue;
            ++flat;
            good += nr.normals.col(i).dot(scan.normals.col(i)) > 0.999;
        }
    }
    CHECK(flat > scan.points.cols() / 2);
    CHECK(good == flat);
}

TEST_CASE("simulator: rays hit the ground and box faces where geometry says", "[sim][ch4]") {
    sim::Scene scene;
    sim::Box b;
    b.center = {10, 0, 1};
    b.half_extent = {1, 2, 1};
    b.id = 7;
    scene.boxes.push_back(b);
    const auto down = sim::cast_ray(scene, {0, 0, 2}, {0, 0, -1}, 0.1, 100);
    CHECK(down.label == -1);
    CHECK(down.range == Approx(2.0));
    const auto ahead = sim::cast_ray(scene, {0, 0, 1}, {1, 0, 0}, 0.1, 100);
    CHECK(ahead.label == 7);
    CHECK(ahead.range == Approx(9.0));
    CHECK(ahead.normal.isApprox(Eigen::Vector3d(-1, 0, 0)));
    // rotate the box by 90 degrees: now its 2 m half-extent faces the ray
    scene.boxes[0].yaw = 1.5707963267948966;
    CHECK(sim::cast_ray(scene, {0, 0, 1}, {1, 0, 0}, 0.1, 100).range == Approx(8.0));
    const auto up = sim::cast_ray(scene, {0, 0, 1}, {0, 0, 1}, 0.1, 100);
    CHECK(up.label == -2);
}

TEST_CASE("PCD: the recorded ASCII scan and the packed binary scans", "[pcd][ch4][data]") {
    const auto a = cloud::read_pcd(std::string(SF_DATA_DIR) + "/parking_lot/run1_scan000_ascii.pcd");
    CHECK(a.fields == std::vector<std::string>{"x", "y", "z", "intensity", "index"});
    CHECK(a.values.cols() == 177);
    CHECK(a.values(0, 0) == Approx(-7.427720069885254));
    CHECK(a.values(4, 0) == 24.0);
    const auto b = cloud::read_pcd(std::string(SF_DATA_DIR) + "/parking_lot/run1_scans.pcd");
    CHECK(b.has("scan"));
    CHECK(b.field("scan").maxCoeff() == 483.0);
    // the first 177 points of the packed file are the ASCII scan, stored as float32
    CHECK((b.xyz().leftCols(177) - a.xyz()).cwiseAbs().maxCoeff() < 1e-6);
}

TEST_CASE("PCD write/read round trip, binary and ASCII", "[pcd][ch4]") {
    cloud::PcdData d;
    d.fields = {"x", "y", "z", "label"};
    d.values.resize(4, 3);
    d.values << 1.5, -2.25, 3.0, 0.1, 0.2, 0.3, -1, 0, 1e3, 0, 1, 2;
    for (bool binary : {true, false}) {
        const std::string path = std::string(binary ? "rt_bin" : "rt_ascii") + ".pcd";
        cloud::write_pcd(path, d, binary, {"label"});
        const auto r = cloud::read_pcd(path);
        CHECK(r.fields == d.fields);
        CHECK((r.values - d.values).cwiseAbs().maxCoeff() < 1e-5);
        std::remove(path.c_str());
    }
}
