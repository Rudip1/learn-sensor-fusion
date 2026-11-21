#pragma once
/// @file fusion.hpp
/// Loosely coupled GNSS + LiDAR fusion: an extended Kalman filter on the planar pose (x, y, yaw) in a local
/// ENU frame, predicted with LiDAR-odometry increments and corrected with GNSS positions, with a chi-square
/// gate on the normalized innovation.
/// Theory: 1_theory/06_gnss_lidar_fusion.md. Equation numbers refer to that file.

#include <Eigen/Core>
#include <vector>

namespace sensor_fusion::fusion {

/// Wrap an angle to (-pi, pi].
double wrap_angle(double a);

/// Pose (x, y, yaw) composed with a body-frame increment (dx, dy, dyaw): eq. (6.1).
Eigen::Vector3d compose(const Eigen::Vector3d& pose, const Eigen::Vector3d& delta);

/// Body-frame increment that takes pose a to pose b (the inverse of compose).
Eigen::Vector3d between(const Eigen::Vector3d& a, const Eigen::Vector3d& b);

/// Covariance of an odometry increment: standard deviations proportional to the distance and the turn,
/// sigma_xy = k_trans |d| + floor_trans, sigma_yaw = k_rot |dyaw| + k_rot_per_m |d| + floor_rot. Eq. (6.3).
struct OdometryNoise {
    double k_trans = 0.02;      ///< translation error per metre travelled
    double k_rot = 0.05;        ///< heading error per radian turned
    double k_rot_per_m = 0.002; ///< heading error per metre travelled [rad/m]
    double floor_trans = 1e-3;  ///< [m]
    double floor_rot = 1e-4;    ///< [rad]
    Eigen::Matrix3d covariance(const Eigen::Vector3d& delta) const;
};

/// EKF on (x, y, yaw).
class PoseEkf {
public:
    PoseEkf(const Eigen::Vector3d& x0, const Eigen::Matrix3d& P0) : x_(x0), P_(P0) {}

    /// Prediction with a body-frame increment and its covariance Q: x <- f(x, u), P <- F P F^T + G Q G^T.
    /// Eqs. (6.1)-(6.2).
    void predict(const Eigen::Vector3d& delta, const Eigen::Matrix3d& Q);

    struct Update {
        Eigen::Vector2d innovation = Eigen::Vector2d::Zero();
        double nis = 0.0;     ///< normalized innovation squared, eq. (6.5)
        bool accepted = false;
    };

    /// Correction with a position measurement z (ENU) of covariance R: eqs. (6.4)-(6.6). The update is skipped
    /// when the NIS exceeds `gate` (pass infinity to accept everything).
    Update update_position(const Eigen::Vector2d& z, const Eigen::Matrix2d& R, double gate);

    const Eigen::Vector3d& state() const { return x_; }
    const Eigen::Matrix3d& covariance() const { return P_; }

private:
    Eigen::Vector3d x_;
    Eigen::Matrix3d P_;
};

/// Chi-square quantile for 2 degrees of freedom: -2 ln(1 - p).
double chi2_2dof(double probability);

struct FusionInput {
    std::vector<double> odom_times;                    ///< time of every odometry pose [s]
    std::vector<Eigen::Vector3d> odom_poses;           ///< odometry poses (x, y, yaw) in the odometry frame
    std::vector<double> gnss_times;                    ///< [s]
    std::vector<Eigen::Vector2d> gnss_positions;       ///< ENU east, north [m]
    std::vector<Eigen::Matrix2d> gnss_covariances;
};

struct FusionOptions {
    OdometryNoise odometry;
    double gate_probability = 0.99;  ///< 0 or >= 1: no gate
    double gnss_scale = 1.0;         ///< multiply the GNSS covariances by this (tuning, section 6.4)
};

struct FusionResult {
    std::vector<double> times;                 ///< one entry per odometry pose
    std::vector<Eigen::Vector3d> states;       ///< filtered (x, y, yaw) after the updates at that time
    std::vector<Eigen::Matrix3d> covariances;
    std::vector<double> nis;                   ///< one per GNSS fix (processed in time order)
    std::vector<bool> accepted;
};

/// Run the filter (Algorithm 6.1): odometry increments between consecutive odometry poses drive the
/// prediction, each GNSS fix updates the state at the latest odometry time not after it.
FusionResult fuse(const FusionInput& input, const Eigen::Vector3d& x0, const Eigen::Matrix3d& P0,
                  const FusionOptions& options = {});

}  // namespace sensor_fusion::fusion
