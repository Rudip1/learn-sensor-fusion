#include "sensor_fusion/registration.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace sensor_fusion::reg {

template <int Dim>
Rigid<Dim> best_rigid(const cloud::Points<Dim>& src, const cloud::Points<Dim>& dst) {
    if (src.cols() != dst.cols() || src.cols() < Dim) throw std::invalid_argument("best_rigid: need >= Dim pairs");
    const cloud::Point<Dim> ms = src.rowwise().mean(), md = dst.rowwise().mean();
    const Eigen::Matrix<double, Dim, Dim> H = (dst.colwise() - md) * (src.colwise() - ms).transpose();  // (5.3)
    const Eigen::JacobiSVD<Eigen::Matrix<double, Dim, Dim>> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix<double, Dim, Dim> S = Eigen::Matrix<double, Dim, Dim>::Identity();
    if ((svd.matrixU() * svd.matrixV().transpose()).determinant() < 0) S(Dim - 1, Dim - 1) = -1;  // no reflection
    Rigid<Dim> T = Rigid<Dim>::Identity();
    T.linear() = svd.matrixU() * S * svd.matrixV().transpose();
    T.translation() = md - T.linear() * ms;
    return T;
}

namespace {

/// Rotation of a small-angle vector: 2-D angle (size 1) or 3-D rotation vector (size 3).
template <int Dim>
Eigen::Matrix<double, Dim, Dim> small_rotation(const Eigen::Matrix<double, Dim == 2 ? 1 : 3, 1>& w);

template <>
Eigen::Matrix2d small_rotation<2>(const Eigen::Matrix<double, 1, 1>& w) {
    return Eigen::Rotation2Dd(w(0)).toRotationMatrix();
}
template <>
Eigen::Matrix3d small_rotation<3>(const Eigen::Vector3d& w) {
    const double a = w.norm();
    return a < 1e-15 ? Eigen::Matrix3d::Identity() : Eigen::AngleAxisd(a, w / a).toRotationMatrix();
}

/// Row of the point-to-plane Jacobian for residual n.(p + w x p + t - q): [d/dw, n^T]. Eq. (5.5).
template <int Dim>
Eigen::Matrix<double, 1, Dim == 2 ? 3 : 6> plane_row(const cloud::Point<Dim>& p, const cloud::Point<Dim>& n);

template <>
Eigen::Matrix<double, 1, 3> plane_row<2>(const Eigen::Vector2d& p, const Eigen::Vector2d& n) {
    Eigen::Matrix<double, 1, 3> J;
    J << n.dot(Eigen::Vector2d(-p.y(), p.x())), n.x(), n.y();  // d(R p)/d theta = perp(p)
    return J;
}
template <>
Eigen::Matrix<double, 1, 6> plane_row<3>(const Eigen::Vector3d& p, const Eigen::Vector3d& n) {
    Eigen::Matrix<double, 1, 6> J;
    J << p.cross(n).transpose(), n.transpose();  // n.(w x p) = w.(p x n)
    return J;
}

}  // namespace

template <int Dim>
IcpResult<Dim> icp(const cloud::Points<Dim>& source, const cloud::Points<Dim>& target, const Rigid<Dim>& initial,
                   const IcpOptions& o) {
    constexpr int R = Dim == 2 ? 1 : 3;  // rotation parameters
    constexpr int P = R + Dim;           // all parameters
    IcpResult<Dim> res;
    res.T = initial;
    if (source.cols() == 0 || target.cols() < Dim + 1) return res;
    const cloud::KdTree<Dim> tree(target);
    cloud::Points<Dim> normals;
    if (o.metric == Metric::PointToPlane) normals = cloud::estimate_normals<Dim>(target, tree, o.normal_k).normals;

    std::vector<int> src_idx, tgt_idx;
    std::vector<double> dist;
    for (int it = 0; it < o.max_iterations; ++it) {
        // 1. correspondences: nearest target point of every transformed source point, gated and trimmed
        const cloud::Points<Dim> moved = cloud::transform<Dim>(res.T, source);
        std::vector<std::pair<double, int>> cand;  // (distance, source index)
        std::vector<int> nn(static_cast<std::size_t>(moved.cols()));
        for (Eigen::Index i = 0; i < moved.cols(); ++i) {
            const auto [j, d2] = tree.nearest(moved.col(i));
            nn[static_cast<std::size_t>(i)] = j;
            if (d2 <= o.max_correspondence * o.max_correspondence) cand.emplace_back(std::sqrt(d2), int(i));
        }
        if (static_cast<int>(cand.size()) < P) break;
        if (o.trim < 1.0) {
            const auto keep = std::max<std::size_t>(std::size_t(P), std::size_t(o.trim * double(cand.size())));
            std::nth_element(cand.begin(), cand.begin() + long(keep) - 1, cand.end());
            cand.resize(keep);
        }
        // 2. minimize over the transform increment
        Rigid<Dim> delta = Rigid<Dim>::Identity();
        double sq = 0.0;
        if (o.metric == Metric::PointToPoint) {
            cloud::Points<Dim> a(Dim, Eigen::Index(cand.size())), b(Dim, Eigen::Index(cand.size()));
            for (std::size_t k = 0; k < cand.size(); ++k) {
                a.col(Eigen::Index(k)) = moved.col(cand[k].second);
                b.col(Eigen::Index(k)) = target.col(nn[std::size_t(cand[k].second)]);
            }
            delta = best_rigid<Dim>(a, b);  // (5.2)
            sq = (cloud::transform<Dim>(delta, a) - b).colwise().squaredNorm().sum();
        } else {
            Eigen::Matrix<double, P, P> A = Eigen::Matrix<double, P, P>::Zero();
            Eigen::Matrix<double, P, 1> g = Eigen::Matrix<double, P, 1>::Zero();
            for (const auto& [d, i] : cand) {
                const cloud::Point<Dim> p = moved.col(i);
                const int j = nn[std::size_t(i)];
                const cloud::Point<Dim> n = normals.col(j);
                const double r = n.dot(p - target.col(j));  // (5.4)
                const auto J = plane_row<Dim>(p, n);
                A += J.transpose() * J;
                g += J.transpose() * r;
            }
            const Eigen::Matrix<double, P, 1> x = -A.ldlt().solve(g);  // (5.5)
            delta.linear() = small_rotation<Dim>(x.template head<R>());
            delta.translation() = x.template tail<Dim>();
            for (const auto& [d, i] : cand) {
                const int j = nn[std::size_t(i)];
                const double r = normals.col(j).dot(delta * cloud::Point<Dim>(moved.col(i)) - target.col(j));
                sq += r * r;
            }
        }
        res.T = delta * res.T;
        res.iterations = it + 1;
        res.pairs = static_cast<int>(cand.size());
        res.rmse = std::sqrt(sq / double(cand.size()));
        res.history.push_back(res.rmse);
        // 3. converged when the increment is negligible
        double angle = 0.0;
        if constexpr (Dim == 2) {
            angle = std::abs(Eigen::Rotation2Dd(delta.linear()).angle());
        } else {
            angle = Eigen::AngleAxisd(delta.linear()).angle();
        }
        if (delta.translation().norm() < o.tolerance && angle < o.tolerance) {
            res.converged = true;
            break;
        }
    }
    return res;
}

template <int Dim>
OdometryResult<Dim> lidar_odometry(const std::vector<cloud::Points<Dim>>& scans, const OdometryOptions& o) {
    OdometryResult<Dim> out;
    if (scans.empty()) return out;
    auto prep = [&](const cloud::Points<Dim>& s) {
        return o.voxel > 0 ? cloud::voxel_downsample<Dim>(s, o.voxel).centroids : s;
    };
    out.poses.push_back(Rigid<Dim>::Identity());
    cloud::Points<Dim> prev = prep(scans[0]);
    Rigid<Dim> motion = Rigid<Dim>::Identity();  // T_{k-1,k}: scan k in the frame of scan k-1
    for (std::size_t k = 1; k < scans.size(); ++k) {
        const cloud::Points<Dim> cur = prep(scans[k]);
        const Rigid<Dim> guess = o.constant_velocity ? motion : Rigid<Dim>::Identity();
        IcpResult<Dim> r = icp<Dim>(cur, prev, guess, o.icp);
        motion = r.T;
        out.poses.push_back(out.poses.back() * motion);  // (5.6)
        out.steps.push_back(std::move(r));
        prev = cur;
    }
    return out;
}

template Rigid<2> best_rigid<2>(const cloud::Points<2>&, const cloud::Points<2>&);
template Rigid<3> best_rigid<3>(const cloud::Points<3>&, const cloud::Points<3>&);
template IcpResult<2> icp<2>(const cloud::Points<2>&, const cloud::Points<2>&, const Rigid<2>&, const IcpOptions&);
template IcpResult<3> icp<3>(const cloud::Points<3>&, const cloud::Points<3>&, const Rigid<3>&, const IcpOptions&);
template OdometryResult<2> lidar_odometry<2>(const std::vector<cloud::Points<2>>&, const OdometryOptions&);
template OdometryResult<3> lidar_odometry<3>(const std::vector<cloud::Points<3>>&, const OdometryOptions&);

}  // namespace sensor_fusion::reg
