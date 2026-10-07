#pragma once
/// @file cloud.hpp
/// Point clouds: PCD file I/O, a k-d tree for nearest-neighbour queries, voxel-grid downsampling and
/// surface normals by local PCA. Templates over the dimension (2 for planar scans, 3 for 3-D LiDAR).
/// Theory: 1_theory/04_point_clouds.md.

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace sensor_fusion::cloud {

template <int Dim>
using Point = Eigen::Matrix<double, Dim, 1>;
/// A point set: one point per column.
template <int Dim>
using Points = Eigen::Matrix<double, Dim, Eigen::Dynamic>;

// ---------------------------------------------------------------------------------------------------
// PCD files (section 4.2)
// ---------------------------------------------------------------------------------------------------

/// Contents of a PCD file: every field as a row of doubles, in file order.
struct PcdData {
    std::vector<std::string> fields;  ///< e.g. {"x", "y", "z", "intensity"}
    Eigen::MatrixXd values;           ///< fields.size() x number of points
    int width = 0;
    int height = 0;
    Eigen::Matrix<double, 7, 1> viewpoint = (Eigen::Matrix<double, 7, 1>() << 0, 0, 0, 1, 0, 0, 0).finished();

    /// Row of a named field; throws if absent.
    Eigen::VectorXd field(const std::string& name) const;
    bool has(const std::string& name) const;
    /// 3 x n matrix of x, y, z (z = 0 if the file has no z field).
    Eigen::Matrix3Xd xyz() const;
};

/// Read an ASCII or binary PCD file (not binary_compressed). Field types F, I, U of sizes 1, 2, 4, 8 and
/// COUNT 1; type letters are accepted in either case.
PcdData read_pcd(const std::string& path);

/// Write a PCD file; every field is stored as float32 ("F 4") except the ones listed in @p uint_fields,
/// which are stored as uint32.
void write_pcd(const std::string& path, const PcdData& data, bool binary = true,
               const std::vector<std::string>& uint_fields = {});

// ---------------------------------------------------------------------------------------------------
// k-d tree (section 4.4)
// ---------------------------------------------------------------------------------------------------

/// Static k-d tree over a copy of the points. Leaves hold up to `leaf_size` points; inner nodes split at
/// the median of the coordinate with the largest spread (Algorithm 4.2).
template <int Dim>
class KdTree {
public:
    explicit KdTree(const Points<Dim>& points, int leaf_size = 10);

    /// Nearest neighbour: (index, squared distance). Index -1 for an empty tree.
    std::pair<int, double> nearest(const Point<Dim>& query) const;
    /// k nearest neighbours sorted by distance (fewer if the tree is smaller).
    void knn(const Point<Dim>& query, int k, std::vector<int>& indices, std::vector<double>& sq_dists) const;
    /// All points within distance r, sorted by distance.
    void radius(const Point<Dim>& query, double r, std::vector<int>& indices, std::vector<double>& sq_dists) const;

    const Points<Dim>& points() const { return points_; }
    int size() const { return static_cast<int>(points_.cols()); }
    int depth() const { return depth_; }

private:
    struct Node {
        int begin = 0, end = 0;          ///< range in index_ (leaves)
        int left = -1, right = -1;       ///< children (inner nodes)
        int axis = -1;                   ///< split coordinate, -1 for a leaf
        double split = 0.0;
    };
    int build(int begin, int end, int depth);

    Points<Dim> points_;
    std::vector<int> index_;
    std::vector<Node> nodes_;
    int leaf_size_;
    int depth_ = 0;
};

/// Brute-force k nearest neighbours, the O(n) reference the k-d tree is tested against.
template <int Dim>
void brute_force_knn(const Points<Dim>& points, const Point<Dim>& query, int k, std::vector<int>& indices,
                     std::vector<double>& sq_dists);

// ---------------------------------------------------------------------------------------------------
// Voxel grid (section 4.3)
// ---------------------------------------------------------------------------------------------------

/// Result of voxel downsampling.
template <int Dim>
struct VoxelResult {
    Points<Dim> centroids;                      ///< one per occupied voxel, in order of first occurrence
    std::vector<int> voxel_of_point;            ///< for every input point, the index of its voxel
    std::vector<int> count;                     ///< points per voxel
    std::vector<Eigen::Matrix<long long, Dim, 1>> keys;  ///< integer voxel coordinates, eq. (4.2)
};

/// Replace all points falling in the same cube of edge @p leaf by their centroid. The grid is anchored at
/// @p origin: voxel key = floor((p - origin) / leaf). Eqs. (4.2)-(4.3).
template <int Dim>
VoxelResult<Dim> voxel_downsample(const Points<Dim>& points, double leaf,
                                  const Point<Dim>& origin = Point<Dim>::Zero());

// ---------------------------------------------------------------------------------------------------
// Normals (section 4.5)
// ---------------------------------------------------------------------------------------------------

/// Centroid and covariance of a point set. Eq. (4.4).
template <int Dim>
std::pair<Point<Dim>, Eigen::Matrix<double, Dim, Dim>> mean_and_covariance(const Points<Dim>& points);

/// Normal of the best-fitting hyperplane (line in 2-D, plane in 3-D): eigenvector of the smallest
/// eigenvalue of the covariance; returns (normal, centroid, surface variation lambda_min / sum lambda).
/// Eqs. (4.4)-(4.6).
template <int Dim>
struct PlaneFit {
    Point<Dim> normal;
    Point<Dim> centroid;
    double variation = 0.0;
};
template <int Dim>
PlaneFit<Dim> fit_plane(const Points<Dim>& points);

template <int Dim>
struct NormalsResult {
    Points<Dim> normals;        ///< unit normals, oriented towards the viewpoint
    Eigen::VectorXd variation;  ///< surface variation per point, 0 on a perfect plane
};

/// Normals of every point from its k nearest neighbours (the point itself included), oriented so that
/// n . (viewpoint - p) >= 0. Algorithm 4.3.
template <int Dim>
NormalsResult<Dim> estimate_normals(const Points<Dim>& points, const KdTree<Dim>& tree, int k,
                                    const Point<Dim>& viewpoint = Point<Dim>::Zero());

// ---------------------------------------------------------------------------------------------------
// Rigid transforms of point sets
// ---------------------------------------------------------------------------------------------------

/// p' = R p + t for every column.
template <int Dim>
Points<Dim> transform(const Eigen::Transform<double, Dim, Eigen::Isometry>& T, const Points<Dim>& points) {
    return (T.linear() * points).colwise() + T.translation();
}

}  // namespace sensor_fusion::cloud
