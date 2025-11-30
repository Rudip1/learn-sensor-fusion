# src/trajectory/fusion_gnss_icp.py
import numpy as np
from .odom_from_gnss import GNSSOdometry

# from .icp_relative_pose import compute_icp_trajectory
from src.trajectory.icp_relative_pose import compute_icp_trajectory


class GNSSICPFusion:
    """
    Fuse GNSS absolute positions with ICP relative/global poses
    to get a globally consistent trajectory with local scan alignment.
    """

    def __init__(self, gnss_loader=None):
        self.odom = GNSSOdometry(loader=gnss_loader)

    def fuse(self, pcd_folder, max_scans=None):
        # 1. Get GNSS poses
        gnss_poses, gnss_positions, _ = self.odom.compute_poses()
        num_scans = len(gnss_poses)
        if max_scans:
            num_scans = min(max_scans, num_scans)

        # 2. Get ICP global poses
        filenames, icp_global_poses = compute_icp_trajectory(
            pcd_folder, max_scans=num_scans
        )

        # 3. Fuse ICP relative deltas on top of GNSS
        fused_poses = [gnss_poses[0]]  # start at GNSS origin

        for i in range(1, num_scans):
            # Compute relative ICP delta from i-1 to i
            icp_delta = np.linalg.inv(icp_global_poses[i - 1]) @ icp_global_poses[i]

            # Apply ICP delta to previous fused pose
            fused_pose = fused_poses[-1] @ icp_delta

            # Optionally correct slightly toward GNSS absolute position
            gnss_pos = gnss_positions[i].flatten()
            fused_pose[:3, 3] = 0.9 * fused_pose[:3, 3] + 0.1 * gnss_pos

            fused_poses.append(fused_pose)

        return fused_poses
