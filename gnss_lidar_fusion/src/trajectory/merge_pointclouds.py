# src/trajectory/merge_pointclouds.py
import numpy as np
import open3d as o3d
from ..utils.transforms import apply_transformation
from ..utils.pcd_loader import load_all_pcds


def merge_pointclouds(pcd_folder, poses):
    """Transform & merge all PCD files into one global point cloud."""
    pcd_dict = load_all_pcds(pcd_folder)
    all_points = []

    for fname, T in zip(sorted(pcd_dict.keys()), poses):
        transformed = apply_transformation(pcd_dict[fname], T)
        all_points.append(transformed)

    return np.vstack(all_points)


def save_xyz(points, filename):
    """Save merged point cloud in plain ASCII XYZ format."""
    np.savetxt(filename, points, fmt="%.6f %.6f %.6f")
    print(f"[OK] Saved XYZ point cloud to: {filename}")
