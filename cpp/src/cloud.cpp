#include "sensor_fusion/cloud.hpp"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace sensor_fusion::cloud {

// ===================================================================================================
// PCD
// ===================================================================================================

Eigen::VectorXd PcdData::field(const std::string& name) const {
    for (std::size_t i = 0; i < fields.size(); ++i)
        if (fields[i] == name) return values.row(static_cast<Eigen::Index>(i)).transpose();
    throw std::out_of_range("PCD has no field '" + name + "'");
}

bool PcdData::has(const std::string& name) const {
    return std::find(fields.begin(), fields.end(), name) != fields.end();
}

Eigen::Matrix3Xd PcdData::xyz() const {
    Eigen::Matrix3Xd p = Eigen::Matrix3Xd::Zero(3, values.cols());
    p.row(0) = field("x").transpose();
    p.row(1) = field("y").transpose();
    if (has("z")) p.row(2) = field("z").transpose();
    return p;
}

namespace {

struct FieldSpec {
    std::string name;
    int size = 4;
    char type = 'F';
    int count = 1;
};

double read_binary_value(const char* p, const FieldSpec& f) {
    switch (f.type) {
        case 'F':
            if (f.size == 4) {
                float v;
                std::memcpy(&v, p, 4);
                return v;
            }
            if (f.size == 8) {
                double v;
                std::memcpy(&v, p, 8);
                return v;
            }
            break;
        case 'I':
            if (f.size == 1) return static_cast<double>(*reinterpret_cast<const std::int8_t*>(p));
            if (f.size == 2) {
                std::int16_t v;
                std::memcpy(&v, p, 2);
                return v;
            }
            if (f.size == 4) {
                std::int32_t v;
                std::memcpy(&v, p, 4);
                return v;
            }
            if (f.size == 8) {
                std::int64_t v;
                std::memcpy(&v, p, 8);
                return static_cast<double>(v);
            }
            break;
        case 'U':
            if (f.size == 1) return static_cast<double>(*reinterpret_cast<const std::uint8_t*>(p));
            if (f.size == 2) {
                std::uint16_t v;
                std::memcpy(&v, p, 2);
                return v;
            }
            if (f.size == 4) {
                std::uint32_t v;
                std::memcpy(&v, p, 4);
                return v;
            }
            if (f.size == 8) {
                std::uint64_t v;
                std::memcpy(&v, p, 8);
                return static_cast<double>(v);
            }
            break;
    }
    throw std::runtime_error("unsupported PCD field type/size for '" + f.name + "'");
}

}  // namespace

PcdData read_pcd(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open PCD file: " + path);

    std::vector<FieldSpec> specs;
    PcdData out;
    long long points = -1;
    std::string data_mode, line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string key;
        ss >> key;
        std::vector<std::string> tok;
        for (std::string t; ss >> t;) tok.push_back(t);
        if (key == "VERSION") continue;
        if (key == "FIELDS") {
            specs.resize(tok.size());
            for (std::size_t i = 0; i < tok.size(); ++i) specs[i].name = tok[i];
        } else if (key == "SIZE") {
            for (std::size_t i = 0; i < tok.size() && i < specs.size(); ++i) specs[i].size = std::stoi(tok[i]);
        } else if (key == "TYPE") {
            for (std::size_t i = 0; i < tok.size() && i < specs.size(); ++i)
                specs[i].type = static_cast<char>(std::toupper(static_cast<unsigned char>(tok[i][0])));
        } else if (key == "COUNT") {
            for (std::size_t i = 0; i < tok.size() && i < specs.size(); ++i) specs[i].count = std::stoi(tok[i]);
        } else if (key == "WIDTH") {
            out.width = std::stoi(tok.at(0));
        } else if (key == "HEIGHT") {
            out.height = std::stoi(tok.at(0));
        } else if (key == "VIEWPOINT") {
            for (std::size_t i = 0; i < 7 && i < tok.size(); ++i) out.viewpoint(static_cast<Eigen::Index>(i)) = std::stod(tok[i]);
        } else if (key == "POINTS") {
            points = std::stoll(tok.at(0));
        } else if (key == "DATA") {
            data_mode = tok.at(0);
            break;
        } else {
            throw std::runtime_error("unexpected PCD header line: " + line);
        }
    }
    if (specs.empty() || data_mode.empty()) throw std::runtime_error("incomplete PCD header: " + path);
    for (const auto& s : specs)
        if (s.count != 1) throw std::runtime_error("PCD fields with COUNT > 1 are not supported ('" + s.name + "')");
    if (points < 0) points = static_cast<long long>(out.width) * std::max(out.height, 1);

    out.fields.clear();
    for (const auto& s : specs) out.fields.push_back(s.name);
    out.values.resize(static_cast<Eigen::Index>(specs.size()), points);

    if (data_mode == "ascii") {
        for (long long i = 0; i < points; ++i)
            for (std::size_t f = 0; f < specs.size(); ++f) {
                std::string tok;
                if (!(in >> tok)) throw std::runtime_error("PCD ascii data ends early: " + path);
                out.values(static_cast<Eigen::Index>(f), i) = (tok == "nan" || tok == "NaN")
                                                                  ? std::numeric_limits<double>::quiet_NaN()
                                                                  : std::stod(tok);
            }
    } else if (data_mode == "binary") {
        int record = 0;
        for (const auto& s : specs) record += s.size;
        std::vector<char> buf(static_cast<std::size_t>(record) * static_cast<std::size_t>(points));
        in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
        if (in.gcount() != static_cast<std::streamsize>(buf.size()))
            throw std::runtime_error("PCD binary data ends early: " + path);
        for (long long i = 0; i < points; ++i) {
            const char* rec = buf.data() + i * record;
            int off = 0;
            for (std::size_t f = 0; f < specs.size(); ++f) {
                out.values(static_cast<Eigen::Index>(f), i) = read_binary_value(rec + off, specs[f]);
                off += specs[f].size;
            }
        }
    } else {
        throw std::runtime_error("PCD DATA mode '" + data_mode + "' is not supported");
    }
    return out;
}

void write_pcd(const std::string& path, const PcdData& data, bool binary,
               const std::vector<std::string>& uint_fields) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write PCD file: " + path);
    const Eigen::Index n = data.values.cols();
    std::vector<bool> is_uint;
    for (const auto& f : data.fields)
        is_uint.push_back(std::find(uint_fields.begin(), uint_fields.end(), f) != uint_fields.end());

    out << "# .PCD v0.7 - Point Cloud Data file format\nVERSION 0.7\nFIELDS";
    for (const auto& f : data.fields) out << ' ' << f;
    out << "\nSIZE";
    for (std::size_t i = 0; i < data.fields.size(); ++i) out << " 4";
    out << "\nTYPE";
    for (std::size_t i = 0; i < data.fields.size(); ++i) out << (is_uint[i] ? " U" : " F");
    out << "\nCOUNT";
    for (std::size_t i = 0; i < data.fields.size(); ++i) out << " 1";
    out << "\nWIDTH " << n << "\nHEIGHT 1\nVIEWPOINT";
    for (int i = 0; i < 7; ++i) out << ' ' << data.viewpoint(i);
    out << "\nPOINTS " << n << "\nDATA " << (binary ? "binary" : "ascii") << "\n";
    for (Eigen::Index i = 0; i < n; ++i) {
        for (std::size_t f = 0; f < data.fields.size(); ++f) {
            const double v = data.values(static_cast<Eigen::Index>(f), i);
            if (binary) {
                if (is_uint[f]) {
                    const auto u = static_cast<std::uint32_t>(v);
                    out.write(reinterpret_cast<const char*>(&u), 4);
                } else {
                    const auto x = static_cast<float>(v);
                    out.write(reinterpret_cast<const char*>(&x), 4);
                }
            } else {
                if (f) out << ' ';
                if (is_uint[f])
                    out << static_cast<std::uint32_t>(v);
                else
                    out << static_cast<float>(v);
            }
        }
        if (!binary) out << '\n';
    }
}

// ===================================================================================================
// k-d tree
// ===================================================================================================

template <int Dim>
KdTree<Dim>::KdTree(const Points<Dim>& points, int leaf_size) : points_(points), leaf_size_(std::max(1, leaf_size)) {
    index_.resize(static_cast<std::size_t>(points_.cols()));
    std::iota(index_.begin(), index_.end(), 0);
    if (!index_.empty()) build(0, static_cast<int>(index_.size()), 1);
}

template <int Dim>
int KdTree<Dim>::build(int begin, int end, int depth) {
    depth_ = std::max(depth_, depth);
    const int id = static_cast<int>(nodes_.size());
    nodes_.push_back(Node{begin, end, -1, -1, -1, 0.0});
    if (end - begin <= leaf_size_) return id;

    // Split along the coordinate with the largest spread, at the median (Algorithm 4.2).
    Point<Dim> lo = points_.col(index_[static_cast<std::size_t>(begin)]), hi = lo;
    for (int i = begin; i < end; ++i) {
        lo = lo.cwiseMin(points_.col(index_[static_cast<std::size_t>(i)]));
        hi = hi.cwiseMax(points_.col(index_[static_cast<std::size_t>(i)]));
    }
    int axis = 0;
    (hi - lo).maxCoeff(&axis);
    if (hi(axis) - lo(axis) <= 0.0) return id;  // all points identical: keep as a leaf

    const int mid = begin + (end - begin) / 2;
    std::nth_element(index_.begin() + begin, index_.begin() + mid, index_.begin() + end,
                     [&](int a, int b) { return points_(axis, a) < points_(axis, b); });
    const double split = points_(axis, index_[static_cast<std::size_t>(mid)]);

    nodes_[static_cast<std::size_t>(id)].axis = axis;
    nodes_[static_cast<std::size_t>(id)].split = split;
    const int left = build(begin, mid, depth + 1);
    const int right = build(mid, end, depth + 1);
    nodes_[static_cast<std::size_t>(id)].left = left;
    nodes_[static_cast<std::size_t>(id)].right = right;
    return id;
}

namespace {

/// Bounded max-heap of (squared distance, index) used by the k-NN search.
struct KBest {
    explicit KBest(int k) : k(k) {}
    double worst() const {
        return static_cast<int>(heap.size()) < k ? std::numeric_limits<double>::infinity() : heap.top().first;
    }
    void offer(double d2, int i) {
        if (static_cast<int>(heap.size()) < k)
            heap.emplace(d2, i);
        else if (d2 < heap.top().first) {
            heap.pop();
            heap.emplace(d2, i);
        }
    }
    int k;
    std::priority_queue<std::pair<double, int>> heap;
};

}  // namespace

template <int Dim>
void KdTree<Dim>::knn(const Point<Dim>& q, int k, std::vector<int>& indices, std::vector<double>& sq_dists) const {
    indices.clear();
    sq_dists.clear();
    if (nodes_.empty() || k <= 0) return;
    KBest best(k);
    // Iterative depth-first search; the far child is visited only if the splitting plane is closer than the
    // current k-th distance.
    std::vector<std::pair<int, double>> stack{{0, 0.0}};  // node, squared distance to its region (lower bound)
    while (!stack.empty()) {
        const auto [id, bound] = stack.back();
        stack.pop_back();
        if (bound >= best.worst()) continue;
        const Node& n = nodes_[static_cast<std::size_t>(id)];
        if (n.axis < 0) {
            for (int i = n.begin; i < n.end; ++i) {
                const int p = index_[static_cast<std::size_t>(i)];
                best.offer((points_.col(p) - q).squaredNorm(), p);
            }
            continue;
        }
        const double diff = q(n.axis) - n.split;
        const int near = diff < 0 ? n.left : n.right;
        const int far = diff < 0 ? n.right : n.left;
        stack.emplace_back(far, std::max(bound, diff * diff));  // pushed first, visited last
        stack.emplace_back(near, bound);
    }
    indices.resize(best.heap.size());
    sq_dists.resize(best.heap.size());
    for (auto i = static_cast<int>(best.heap.size()) - 1; i >= 0; --i) {
        sq_dists[static_cast<std::size_t>(i)] = best.heap.top().first;
        indices[static_cast<std::size_t>(i)] = best.heap.top().second;
        best.heap.pop();
    }
}

template <int Dim>
std::pair<int, double> KdTree<Dim>::nearest(const Point<Dim>& q) const {
    std::vector<int> idx;
    std::vector<double> d2;
    knn(q, 1, idx, d2);
    if (idx.empty()) return {-1, std::numeric_limits<double>::infinity()};
    return {idx[0], d2[0]};
}

template <int Dim>
void KdTree<Dim>::radius(const Point<Dim>& q, double r, std::vector<int>& indices,
                         std::vector<double>& sq_dists) const {
    indices.clear();
    sq_dists.clear();
    if (nodes_.empty()) return;
    const double r2 = r * r;
    std::vector<std::pair<double, int>> found;
    std::vector<int> stack{0};
    while (!stack.empty()) {
        const Node& n = nodes_[static_cast<std::size_t>(stack.back())];
        stack.pop_back();
        if (n.axis < 0) {
            for (int i = n.begin; i < n.end; ++i) {
                const int p = index_[static_cast<std::size_t>(i)];
                const double d2 = (points_.col(p) - q).squaredNorm();
                if (d2 <= r2) found.emplace_back(d2, p);
            }
            continue;
        }
        const double diff = q(n.axis) - n.split;
        if (diff - r <= 0) stack.push_back(n.left);  // ball reaches the lower side
        if (diff + r >= 0) stack.push_back(n.right);  // ball reaches the upper side
    }
    std::sort(found.begin(), found.end());
    for (const auto& [d2, p] : found) {
        sq_dists.push_back(d2);
        indices.push_back(p);
    }
}

template <int Dim>
void brute_force_knn(const Points<Dim>& points, const Point<Dim>& q, int k, std::vector<int>& indices,
                     std::vector<double>& sq_dists) {
    std::vector<std::pair<double, int>> all(static_cast<std::size_t>(points.cols()));
    for (Eigen::Index i = 0; i < points.cols(); ++i)
        all[static_cast<std::size_t>(i)] = {(points.col(i) - q).squaredNorm(), static_cast<int>(i)};
    const auto kk = std::min<std::size_t>(static_cast<std::size_t>(std::max(k, 0)), all.size());
    std::partial_sort(all.begin(), all.begin() + static_cast<std::ptrdiff_t>(kk), all.end());
    indices.resize(kk);
    sq_dists.resize(kk);
    for (std::size_t i = 0; i < kk; ++i) {
        sq_dists[i] = all[i].first;
        indices[i] = all[i].second;
    }
}

// ===================================================================================================
// Voxel grid
// ===================================================================================================

namespace {

template <int Dim>
struct KeyHash {
    std::size_t operator()(const Eigen::Matrix<long long, Dim, 1>& k) const {
        std::size_t h = 1469598103934665603ULL;
        for (int i = 0; i < Dim; ++i) {
            h ^= static_cast<std::size_t>(k(i)) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return h;
    }
};

}  // namespace

template <int Dim>
VoxelResult<Dim> voxel_downsample(const Points<Dim>& points, double leaf, const Point<Dim>& origin) {
    if (!(leaf > 0.0)) throw std::invalid_argument("voxel_downsample: leaf size must be positive");
    using Key = Eigen::Matrix<long long, Dim, 1>;
    VoxelResult<Dim> r;
    std::unordered_map<Key, int, KeyHash<Dim>> lookup;
    std::vector<Point<Dim>> sums;
    r.voxel_of_point.resize(static_cast<std::size_t>(points.cols()));
    for (Eigen::Index i = 0; i < points.cols(); ++i) {
        const Key key = ((points.col(i) - origin) / leaf).array().floor().template cast<long long>();  // eq. (4.2)
        auto [it, inserted] = lookup.try_emplace(key, static_cast<int>(sums.size()));
        if (inserted) {
            sums.push_back(Point<Dim>::Zero());
            r.count.push_back(0);
            r.keys.push_back(key);
        }
        sums[static_cast<std::size_t>(it->second)] += points.col(i);
        ++r.count[static_cast<std::size_t>(it->second)];
        r.voxel_of_point[static_cast<std::size_t>(i)] = it->second;
    }
    r.centroids.resize(Dim, static_cast<Eigen::Index>(sums.size()));
    for (std::size_t v = 0; v < sums.size(); ++v)
        r.centroids.col(static_cast<Eigen::Index>(v)) = sums[v] / r.count[v];  // eq. (4.3)
    return r;
}

// ===================================================================================================
// Normals
// ===================================================================================================

template <int Dim>
std::pair<Point<Dim>, Eigen::Matrix<double, Dim, Dim>> mean_and_covariance(const Points<Dim>& points) {
    const Point<Dim> mean = points.rowwise().mean();
    const Points<Dim> d = points.colwise() - mean;
    const Eigen::Matrix<double, Dim, Dim> cov = d * d.transpose() / static_cast<double>(std::max<Eigen::Index>(1, points.cols()));
    return {mean, cov};  // eq. (4.4)
}

template <int Dim>
PlaneFit<Dim> fit_plane(const Points<Dim>& points) {
    const auto [mean, cov] = mean_and_covariance<Dim>(points);
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, Dim, Dim>> es(cov);
    PlaneFit<Dim> f;
    f.centroid = mean;
    f.normal = es.eigenvectors().col(0);  // smallest eigenvalue first, eq. (4.5)
    const double sum = es.eigenvalues().sum();
    f.variation = sum > 0.0 ? es.eigenvalues()(0) / sum : 0.0;  // eq. (4.6)
    return f;
}

template <int Dim>
NormalsResult<Dim> estimate_normals(const Points<Dim>& points, const KdTree<Dim>& tree, int k,
                                    const Point<Dim>& viewpoint) {
    NormalsResult<Dim> r;
    r.normals.resize(Dim, points.cols());
    r.variation.resize(points.cols());
    std::vector<int> idx;
    std::vector<double> d2;
    Points<Dim> nb;
    for (Eigen::Index i = 0; i < points.cols(); ++i) {
        tree.knn(points.col(i), k, idx, d2);
        nb.resize(Dim, static_cast<Eigen::Index>(idx.size()));
        for (std::size_t j = 0; j < idx.size(); ++j) nb.col(static_cast<Eigen::Index>(j)) = tree.points().col(idx[j]);
        PlaneFit<Dim> f = fit_plane<Dim>(nb);
        if (f.normal.dot(viewpoint - points.col(i)) < 0.0) f.normal = -f.normal;  // orient to the sensor
        r.normals.col(i) = f.normal;
        r.variation(i) = f.variation;
    }
    return r;
}

// Explicit instantiations for planar and 3-D clouds.
template class KdTree<2>;
template class KdTree<3>;
template void brute_force_knn<2>(const Points<2>&, const Point<2>&, int, std::vector<int>&, std::vector<double>&);
template void brute_force_knn<3>(const Points<3>&, const Point<3>&, int, std::vector<int>&, std::vector<double>&);
template VoxelResult<2> voxel_downsample<2>(const Points<2>&, double, const Point<2>&);
template VoxelResult<3> voxel_downsample<3>(const Points<3>&, double, const Point<3>&);
template std::pair<Point<2>, Eigen::Matrix2d> mean_and_covariance<2>(const Points<2>&);
template std::pair<Point<3>, Eigen::Matrix3d> mean_and_covariance<3>(const Points<3>&);
template PlaneFit<2> fit_plane<2>(const Points<2>&);
template PlaneFit<3> fit_plane<3>(const Points<3>&);
template NormalsResult<2> estimate_normals<2>(const Points<2>&, const KdTree<2>&, int, const Point<2>&);
template NormalsResult<3> estimate_normals<3>(const Points<3>&, const KdTree<3>&, int, const Point<3>&);

}  // namespace sensor_fusion::cloud
