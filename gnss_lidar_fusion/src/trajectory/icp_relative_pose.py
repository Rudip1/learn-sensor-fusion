import open3d as o3d
import numpy as np
from typing import List, Tuple
from src.utils.pcd_loader import load_all_pcds


def compute_icp_trajectory(
    pcd_folder: str,
    max_scans: int | None = None,
    voxel_size: float = 0.2,
    max_correspondence_distance: float = 1.0,
    gnss_positions: np.ndarray | None = None,  # GNSS guidance
) -> Tuple[List[str], List[np.ndarray]]:
    """
    Computes the global trajectory using ICP between consecutive scans,
    optionally aligned with GNSS starting direction.

    Args:
        pcd_folder: folder containing PCD files.
        max_scans: max number of scans to process.
        voxel_size: voxel size for downsampling point clouds.
        max_correspondence_distance: max distance threshold for ICP.
        gnss_positions: optional GNSS positions to align initial ICP direction.

    Returns:
        filenames: list of PCD filenames used.
        poses: list of 4x4 numpy arrays representing global poses.
    """
    # Load all PCDs
    pcd_dict = load_all_pcds(pcd_folder)
    if max_scans is None:
        max_scans = len(pcd_dict)
    filenames = sorted(pcd_dict.keys())[:max_scans]

    poses = [np.eye(4)]  # start at origin

    # Prepare first scan
    prev_pcd = o3d.geometry.PointCloud()
    prev_pcd.points = o3d.utility.Vector3dVector(pcd_dict[filenames[0]])
    prev_pcd = prev_pcd.voxel_down_sample(voxel_size)

    for i, fname in enumerate(filenames[1:]):
        curr_pcd = o3d.geometry.PointCloud()
        curr_pcd.points = o3d.utility.Vector3dVector(pcd_dict[fname])
        curr_pcd = curr_pcd.voxel_down_sample(voxel_size)

        # Initial guess: use GNSS delta if provided
        if gnss_positions is not None:
            delta = gnss_positions[i + 1] - gnss_positions[i]
            trans_init = np.eye(4)
            trans_init[:3, 3] = delta  # translate along GNSS direction
        else:
            trans_init = np.eye(4)

        # ICP point-to-point
        reg_p2p = o3d.pipelines.registration.registration_icp(
            curr_pcd,
            prev_pcd,
            max_correspondence_distance,
            trans_init,
            o3d.pipelines.registration.TransformationEstimationPointToPoint(),
        )

        # Compute global pose
        global_pose = poses[-1] @ reg_p2p.transformation
        poses.append(global_pose)

        # Update previous point cloud
        prev_pcd = curr_pcd

    return filenames, poses
