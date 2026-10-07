# src/trajectory/odom_from_gnss.py
from typing import Optional
import numpy as np
from src.utils.transforms import create_transformation_matrix
from .gnss_loader import GNSSLoader


class GNSSOdometry:
    """
    Convert GNSS (lat, lon, alt) -> local ENU and build 4x4 pose matrices.
    Uses simple local projection (meters_per_deg approximation).
    """

    def __init__(self, loader: Optional[GNSSLoader] = None):
        self.loader = loader or GNSSLoader()

    @staticmethod
    def latlon_to_local_meters(lat, lon, lat0, lon0):
        """
        Approximate conversion: differences in degrees -> meters.
        x = east, y = north (we'll place north as Y), z = up.
        Uses:
          meters per degree lat ≈ 111132
          meters per degree lon ≈ 111320 * cos(lat0)
        """
        meters_per_deg_lat = 111132.0
        meters_per_deg_lon = 111320.0 * np.cos(np.deg2rad(lat0))
        east = (lon - lon0) * meters_per_deg_lon
        north = (lat - lat0) * meters_per_deg_lat
        return east, north

    def compute_poses(self, csv_path: Optional[str] = None):
        """
        Returns:
          poses: list of 4x4 numpy arrays (poses in local ENU meters)
          positions: Nx3 numpy array of [x_east, y_north, z_up]
          timestamps_ns: 1D numpy array of header_stamp_secs*1e9 + header_stamp_nsecs
        """
        df = self.loader.load(csv_path)
        # require that csv has the columns you provided
        # convert timestamps to ns for matching later
        ts_ns = (
            df["header_stamp_secs"].astype(np.int64) * 1_000_000_000
            + df["header_stamp_nsecs"].astype(np.int64)
        ).to_numpy()

        lat = df["latitude"].astype(float).to_numpy()
        lon = df["longitude"].astype(float).to_numpy()
        alt = df["altitude"].astype(float).to_numpy()

        lat0 = lat[0]
        lon0 = lon[0]
        alt0 = alt[0]

        positions = []
        poses = []
        for la, lo, al in zip(lat, lon, alt):
            east, north = self.latlon_to_local_meters(la, lo, lat0, lon0)
            up = al - alt0
            # choose ENU ordering for translation: [east, north, up]
            t = np.array([east, north, up]).reshape(3, 1)
            T = create_transformation_matrix(np.eye(3), t)
            poses.append(T)
            positions.append([east, north, up])

        return poses, np.array(positions), ts_ns


if __name__ == "__main__":
    print("Testing GNSS Odometry...")
    odom = GNSSOdometry()
    poses, positions, ts = odom.compute_poses()
    print("Loaded poses:", len(poses))
    print("First pose:\n", poses[0])
    print("First position:", positions[0])


"""
python -m src.trajectory.odom_from_gnss

pravin@pravin-nitro:~/gnss_lidar_fusion$ python -m src.trajectory.odom_from_gnss
Testing GNSS Odometry...
Loaded poses: 662
First pose:
 [[1. 0. 0. 0.]
 [0. 1. 0. 0.]
 [0. 0. 1. 0.]
 [0. 0. 0. 1.]]
First position: [0. 0. 0.]
"""
