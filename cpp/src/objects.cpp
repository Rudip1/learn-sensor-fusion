#include "sensor_fusion/objects.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace sensor_fusion::objects {

GroundResult segment_ground(const Eigen::Matrix3Xd& P, double threshold, int iterations, double max_tilt,
                            std::uint32_t seed) {
    GroundResult r;
    const Eigen::Index n = P.cols();
    r.ground.assign(std::size_t(n), false);
    if (n < 3) return r;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<Eigen::Index> pick(0, n - 1);
    auto count = [&](const Plane& pl) {
        int c = 0;
        for (Eigen::Index i = 0; i < n; ++i) c += std::abs(pl.distance(P.col(i))) < threshold;  // (8.1)
        return c;
    };
    int best = -1;
    for (int it = 0; it < iterations; ++it) {
        const Eigen::Vector3d a = P.col(pick(gen)), b = P.col(pick(gen)), c = P.col(pick(gen));
        Eigen::Vector3d nrm = (b - a).cross(c - a);
        if (nrm.norm() < 1e-9) continue;
        nrm.normalize();
        if (nrm.z() < 0) nrm = -nrm;
        if (std::acos(std::min(1.0, nrm.z())) > max_tilt) continue;  // walls are not ground
        const Plane pl{nrm, -nrm.dot(a)};
        const int k = count(pl);
        if (k > best) best = k, r.plane = pl;
    }
    if (best < 3) return r;
    // refit by PCA to the inliers (chapter 4), then classify
    std::vector<int> in;
    for (Eigen::Index i = 0; i < n; ++i)
        if (std::abs(r.plane.distance(P.col(i))) < threshold) in.push_back(int(i));
    cloud::Points<3> Q(3, Eigen::Index(in.size()));
    for (std::size_t k = 0; k < in.size(); ++k) Q.col(Eigen::Index(k)) = P.col(in[k]);
    const auto fit = cloud::fit_plane<3>(Q);
    Eigen::Vector3d nrm = fit.normal.z() < 0 ? Eigen::Vector3d(-fit.normal) : Eigen::Vector3d(fit.normal);
    r.plane = Plane{nrm, -nrm.dot(fit.centroid)};
    for (Eigen::Index i = 0; i < n; ++i) {
        r.ground[std::size_t(i)] = std::abs(r.plane.distance(P.col(i))) < threshold;
        r.inliers += r.ground[std::size_t(i)];
    }
    return r;
}

std::vector<bool> ground_by_height(const Eigen::Matrix3Xd& P, double z_max) {
    std::vector<bool> g(std::size_t(P.cols()));
    for (Eigen::Index i = 0; i < P.cols(); ++i) g[std::size_t(i)] = P(2, i) < z_max;
    return g;
}

template <int Dim>
std::vector<int> dbscan(const cloud::Points<Dim>& P, double eps, int min_points) {
    const int n = int(P.cols());
    std::vector<int> label(std::size_t(n), -2);  // -2 unvisited, -1 noise
    if (n == 0) return {};
    const cloud::KdTree<Dim> tree(P);
    std::vector<int> nb;
    std::vector<double> d2;
    int cluster = 0;
    for (int i = 0; i < n; ++i) {
        if (label[std::size_t(i)] != -2) continue;
        tree.radius(P.col(i), eps, nb, d2);  // (8.2): the eps-neighbourhood, i itself included
        if (int(nb.size()) < min_points) {
            label[std::size_t(i)] = -1;  // noise for now; may become a border point later
            continue;
        }
        label[std::size_t(i)] = cluster;
        std::vector<int> queue(nb.begin(), nb.end());
        for (std::size_t q = 0; q < queue.size(); ++q) {
            const int j = queue[q];
            if (label[std::size_t(j)] == -1) label[std::size_t(j)] = cluster;  // border point
            if (label[std::size_t(j)] != -2) continue;
            label[std::size_t(j)] = cluster;
            tree.radius(P.col(j), eps, nb, d2);
            if (int(nb.size()) >= min_points) queue.insert(queue.end(), nb.begin(), nb.end());  // core: expand
        }
        ++cluster;
    }
    return label;
}

std::vector<Cluster> make_clusters(const Eigen::Matrix3Xd& P, const std::vector<int>& labels) {
    int k = 0;
    for (int l : labels) k = std::max(k, l + 1);
    std::vector<Cluster> out(static_cast<std::size_t>(k));
    for (std::size_t i = 0; i < labels.size(); ++i)
        if (labels[i] >= 0) out[std::size_t(labels[i])].indices.push_back(int(i));
    for (Cluster& c : out) {
        if (c.indices.empty()) continue;
        c.min = c.max = P.col(c.indices[0]);
        for (int i : c.indices) {
            c.centroid += P.col(i);
            c.min = c.min.cwiseMin(P.col(i));
            c.max = c.max.cwiseMax(P.col(i));
        }
        c.centroid /= double(c.indices.size());
    }
    return out;
}

namespace {
struct KeyHash {
    std::size_t operator()(const Eigen::Vector3i& k) const {
        return std::size_t(k.x()) * 73856093u ^ std::size_t(k.y()) * 19349663u ^ std::size_t(k.z()) * 83492791u;
    }
};
struct KeyEq {
    bool operator()(const Eigen::Vector3i& a, const Eigen::Vector3i& b) const { return a == b; }
};
Eigen::Vector3i key(const Eigen::Vector3d& p, double voxel) {
    return (p / voxel).array().floor().cast<int>();
}
}  // namespace

std::vector<Eigen::VectorXd> persistence(const std::vector<Eigen::Matrix3Xd>& scans, double voxel, bool dilate) {
    // number of scans occupying each voxel
    std::unordered_map<Eigen::Vector3i, int, KeyHash, KeyEq> occ;
    for (const auto& s : scans) {
        std::unordered_set<Eigen::Vector3i, KeyHash, KeyEq> seen;
        for (Eigen::Index i = 0; i < s.cols(); ++i) seen.insert(key(s.col(i), voxel));
        for (const auto& k : seen) ++occ[k];
    }
    // a point's persistence: occupancy of its voxel (optionally the most occupied of it and its 26 neighbours)
    std::vector<Eigen::VectorXd> out;
    const double N = double(scans.size());
    for (const auto& s : scans) {
        Eigen::VectorXd p(s.cols());
        for (Eigen::Index i = 0; i < s.cols(); ++i) {
            const Eigen::Vector3i k = key(s.col(i), voxel);
            int m = 0;
            const int r = dilate ? 1 : 0;
            for (int dx = -r; dx <= r; ++dx)
                for (int dy = -r; dy <= r; ++dy)
                    for (int dz = -r; dz <= r; ++dz) {
                        const auto it = occ.find(k + Eigen::Vector3i(dx, dy, dz));
                        if (it != occ.end()) m = std::max(m, it->second);
                    }
            p(i) = m / N;  // (8.3)
        }
        out.push_back(p);
    }
    return out;
}

std::vector<std::vector<bool>> dynamic_points(const std::vector<Eigen::Matrix3Xd>& scans,
                                              const std::vector<std::vector<int>>& labels, double voxel,
                                              double static_ratio, bool dilate) {
    if (scans.size() != labels.size()) throw std::invalid_argument("dynamic_points: one label vector per scan");
    const auto pers = persistence(scans, voxel, dilate);
    std::vector<std::vector<bool>> out;
    for (std::size_t k = 0; k < scans.size(); ++k) {
        int nc = 0;
        for (int l : labels[k]) nc = std::max(nc, l + 1);
        std::vector<double> sum(std::size_t(nc), 0.0);
        std::vector<int> cnt(std::size_t(nc), 0);
        for (std::size_t i = 0; i < labels[k].size(); ++i)
            if (labels[k][i] >= 0) sum[std::size_t(labels[k][i])] += pers[k](Eigen::Index(i)), ++cnt[std::size_t(labels[k][i])];
        std::vector<bool> dyn(labels[k].size(), false);
        for (std::size_t i = 0; i < labels[k].size(); ++i) {
            const int l = labels[k][i];
            if (l >= 0) dyn[i] = sum[std::size_t(l)] / cnt[std::size_t(l)] < static_ratio;  // (8.4)
        }
        out.push_back(dyn);
    }
    return out;
}

template std::vector<int> dbscan<2>(const cloud::Points<2>&, double, int);
template std::vector<int> dbscan<3>(const cloud::Points<3>&, double, int);

}  // namespace sensor_fusion::objects
