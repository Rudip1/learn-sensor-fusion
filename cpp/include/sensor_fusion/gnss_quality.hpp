#pragma once
/// @file gnss_quality.hpp
/// Dilution of precision from satellite geometry, error ellipses and precision statistics, fix-type codes
/// and a standstill / motion classifier for GNSS tracks.
/// Theory: 1_theory/03_gnss_quality.md.

#include <Eigen/Core>
#include <string>
#include <vector>

namespace sensor_fusion::gnss {

/// Dilution-of-precision values: square roots of sums of diagonal entries of Q = (G^T W G)^-1 in ENU.
struct Dop {
    double gdop = 0.0;  ///< sqrt(q_EE + q_NN + q_UU + q_tt)
    double pdop = 0.0;  ///< sqrt(q_EE + q_NN + q_UU)
    double hdop = 0.0;  ///< sqrt(q_EE + q_NN)
    double vdop = 0.0;  ///< sqrt(q_UU)
    double tdop = 0.0;  ///< sqrt(q_tt)
    Eigen::Matrix4d cofactor_enu = Eigen::Matrix4d::Zero();  ///< Q in (east, north, up, clock)
    bool valid = false;  ///< false with fewer than four satellites or a singular geometry
};

/// Geometry matrix in the local ENU frame from azimuth / elevation [rad]: rows
/// [-cos(el) sin(az), -cos(el) cos(az), -sin(el), 1]. Eq. (3.1).
Eigen::MatrixX4d geometry_matrix_enu(const Eigen::Ref<const Eigen::VectorXd>& azimuth,
                                     const Eigen::Ref<const Eigen::VectorXd>& elevation);

/// DOP from azimuth / elevation [rad]; optional per-satellite weights (empty: unit weights). Eqs. (3.2)-(3.3).
Dop compute_dop(const Eigen::Ref<const Eigen::VectorXd>& azimuth, const Eigen::Ref<const Eigen::VectorXd>& elevation,
                const Eigen::Ref<const Eigen::VectorXd>& weights = Eigen::VectorXd());

/// DOP from an ECEF cofactor matrix (as returned by solve_position) by rotating its position block into
/// ENU at (lat, lon) [rad]. Eq. (3.4).
Dop dop_from_ecef_cofactor(const Eigen::Matrix4d& cofactor_ecef, double lat, double lon);

/// Elevation-dependent weight w = sin^2(el), a common simple model (sigma proportional to 1/sin(el)).
Eigen::VectorXd elevation_weights(const Eigen::Ref<const Eigen::VectorXd>& elevation);

/// Confidence ellipse of a 2-D Gaussian with covariance @p cov. Eq. (3.7).
struct ErrorEllipse {
    double semi_major = 0.0;  ///< [m]
    double semi_minor = 0.0;  ///< [m]
    double angle = 0.0;       ///< direction of the major axis, counter-clockwise from the first axis (east) [rad]
};

/// @param probability probability mass inside the ellipse, e.g. 0.3935 (1-sigma), 0.95.
ErrorEllipse error_ellipse(const Eigen::Matrix2d& cov, double probability = 0.3934693402873666);

/// Precision summary of a set of horizontal positions (static receiver). Eq. (3.6).
struct PrecisionStats {
    Eigen::Vector2d mean = Eigen::Vector2d::Zero();
    Eigen::Matrix2d covariance = Eigen::Matrix2d::Zero();  ///< sample covariance (n - 1)
    double sigma_east = 0.0;
    double sigma_north = 0.0;
    double drms = 0.0;      ///< sqrt(sigma_e^2 + sigma_n^2)
    double two_drms = 0.0;  ///< 2 DRMS
    double cep50 = 0.0;     ///< empirical median of the horizontal distance to the mean
    double r95 = 0.0;       ///< empirical 95th percentile of that distance
    int count = 0;
};

/// @param en 2 x n horizontal positions (east, north) [m].
PrecisionStats precision_stats(const Eigen::Ref<const Eigen::Matrix2Xd>& en);

/// Human-readable names of the fix codes found in NMEA (GGA quality, GSA fix type) and ROS NavSatStatus.
std::string gga_quality_name(int quality);
std::string gsa_fix_type_name(int fix_type);
std::string navsat_status_name(int status);

/// Parameters of the motion classifier (Algorithm 3.1).
struct MotionParams {
    double window_s = 5.0;            ///< displacement is measured over this time span
    /// 1-sigma per-axis noise of the *displacement* between two fixes window_s apart, for a receiver that
    /// stands still. For white noise of sigma per fix it is sqrt(2) sigma; real GNSS errors are correlated
    /// in time, so calibrate it from a static log (eq. 3.9).
    double sigma_displacement_m = 1.0;
    double chi2_threshold = 13.8155;  ///< chi^2_2 quantile; 13.8 -> p = 0.999
    double speed_threshold_mps = 0.5; ///< receiver (Doppler) speed above this counts as motion
    int min_run = 3;                  ///< label changes must persist this many epochs (hysteresis)
};

struct MotionResult {
    std::vector<bool> moving;  ///< per-epoch decision after hysteresis
    Eigen::VectorXd test_statistic;  ///< d^2 of eq. (3.8) for each epoch (NaN when no partner epoch)
    Eigen::VectorXd displacement_speed;  ///< |p(t) - p(t - window)| / window [m/s]
};

/// Classify each epoch as moving or standing. @p t seconds, @p en 2 x n positions, @p speed receiver speed
/// (may be empty or contain NaN). Eq. (3.8), Algorithm 3.1.
MotionResult classify_motion(const Eigen::Ref<const Eigen::VectorXd>& t, const Eigen::Ref<const Eigen::Matrix2Xd>& en,
                             const Eigen::Ref<const Eigen::VectorXd>& speed, const MotionParams& params = MotionParams());

}  // namespace sensor_fusion::gnss
