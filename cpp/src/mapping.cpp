#include "sensor_fusion/mapping.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#include "sensor_fusion/cloud.hpp"

namespace sensor_fusion::mapping {

Eigen::Isometry3d sensor_pose(const Eigen::Vector3d& p, const Eigen::Isometry3d& T_vehicle_sensor) {
    Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
    T.linear() = Eigen::AngleAxisd(p(2), Eigen::Vector3d::UnitZ()).toRotationMatrix();
    T.translation() = Eigen::Vector3d(p(0), p(1), 0.0);
    return T * T_vehicle_sensor;  // (7.1)
}

Map aggregate(const std::vector<Eigen::Matrix3Xd>& scans, const std::vector<Eigen::Isometry3d>& poses, double voxel) {
    if (scans.size() != poses.size()) throw std::invalid_argument("aggregate: one pose per scan");
    Eigen::Index n = 0;
    for (const auto& s : scans) n += s.cols();
    Map m;
    m.points.resize(3, n);
    m.scan.resize(n);
    Eigen::Index c = 0;
    for (std::size_t k = 0; k < scans.size(); ++k) {
        m.points.middleCols(c, scans[k].cols()) = (poses[k].linear() * scans[k]).colwise() + poses[k].translation();  // (7.2)
        m.scan.segment(c, scans[k].cols()).setConstant(int(k));
        c += scans[k].cols();
    }
    if (voxel <= 0) return m;
    const auto v = cloud::voxel_downsample<3>(m.points, voxel);
    Map d;
    d.points = v.centroids;
    d.scan = Eigen::VectorXi::Constant(v.centroids.cols(), -1);
    for (Eigen::Index i = 0; i < m.points.cols(); ++i) {
        int& s = d.scan(v.voxel_of_point[std::size_t(i)]);
        if (s < 0) s = m.scan(i);
    }
    return d;
}

double cloud_to_cloud_rms(const Eigen::Matrix3Xd& a, const Eigen::Matrix3Xd& ref) {
    const cloud::KdTree<3> tree(ref);
    double s = 0;
    for (Eigen::Index i = 0; i < a.cols(); ++i) s += tree.nearest(a.col(i)).second;
    return std::sqrt(s / double(a.cols()));
}

Projection project(const Eigen::Matrix3Xd& world, const Eigen::Isometry3d& T_wc, const sim::Camera& cam) {
    Projection p;
    const Eigen::Matrix3Xd Xc = T_wc.inverse() * world;
    const Eigen::Matrix3Xd x = cam.K * Xc;
    p.pixels.resize(2, world.cols());
    p.depth = Xc.row(2).transpose();
    p.inside.resize(std::size_t(world.cols()));
    for (Eigen::Index i = 0; i < world.cols(); ++i) {
        const double z = Xc(2, i);
        p.pixels.col(i) = x.col(i).head<2>() / z;  // (7.3)
        p.inside[std::size_t(i)] = z > 1e-6 && p.pixels(0, i) >= 0 && p.pixels(1, i) >= 0 &&
                                   p.pixels(0, i) < cam.width && p.pixels(1, i) < cam.height;
    }
    return p;
}

Colors colorize(const Eigen::Matrix3Xd& world, const std::vector<sim::RgbImage>& images,
                const std::vector<Eigen::Isometry3d>& T_wc, const sim::Camera& cam, bool occlusion_check,
                double tolerance, int window) {
    if (images.size() != T_wc.size()) throw std::invalid_argument("colorize: one pose per image");
    Colors c;
    c.rgb = Eigen::Matrix3Xd::Zero(3, world.cols());
    c.camera.assign(std::size_t(world.cols()), -1);
    Eigen::VectorXd best = Eigen::VectorXd::Constant(world.cols(), std::numeric_limits<double>::infinity());
    const double inf = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < images.size(); ++k) {
        const Projection p = project(world, T_wc[k], cam);
        // z-buffer of the map: smallest depth per pixel
        Eigen::MatrixXd zbuf = Eigen::MatrixXd::Constant(cam.height, cam.width, inf);
        if (occlusion_check)
            for (Eigen::Index i = 0; i < world.cols(); ++i)
                if (p.inside[std::size_t(i)]) {
                    double& z = zbuf(int(p.pixels(1, i)), int(p.pixels(0, i)));
                    z = std::min(z, p.depth(i));
                }
        for (Eigen::Index i = 0; i < world.cols(); ++i) {
            if (!p.inside[std::size_t(i)] || p.depth(i) >= best(i)) continue;
            const int u = int(p.pixels(0, i)), v = int(p.pixels(1, i));
            if (occlusion_check) {
                double zmin = inf;
                for (int dv = -window; dv <= window; ++dv)
                    for (int du = -window; du <= window; ++du) {
                        const int uu = u + du, vv = v + dv;
                        if (uu >= 0 && vv >= 0 && uu < cam.width && vv < cam.height) zmin = std::min(zmin, zbuf(vv, uu));
                    }
                if (p.depth(i) > zmin + tolerance) continue;  // hidden behind a closer map point
            }
            best(i) = p.depth(i);
            c.rgb.col(i) = images[k].rgb.col(Eigen::Index(v) * cam.width + u).cast<double>();
            c.camera[std::size_t(i)] = int(k);
        }
    }
    return c;
}

}  // namespace sensor_fusion::mapping
