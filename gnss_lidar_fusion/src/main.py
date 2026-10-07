# src/main.py
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Trajectory methods
from src.trajectory.odom_from_gnss import GNSSOdometry
from src.trajectory.icp_relative_pose import compute_icp_trajectory
from src.trajectory.fusion_gnss_icp import GNSSICPFusion

# Merge / save pointcloud helpers
from src.trajectory.merge_pointclouds import merge_pointclouds, save_xyz

# Utils
from src.utils.transforms import apply_transformation
from src.utils.pcd_loader import get_default_pcd_folder, load_all_pcds

# Visualization tools (expects the visualization module you've prepared)
from src.utils.visualization import (
    visualize_trajectory_2d,
    visualize_trajectory_3d,
    visualize_trajectories_3d_comparison,  # GNSS vs Fused 3D comparator
    visualize_trajectories_2d_comparison as vis2d_from_vis,
)


# 2D helper: GNSS vs FUSED
def visualize_trajectories_2d_comparison(
    gnss_positions, fused_positions, title="GNSS vs Fused (Top-down)"
):
    """2D comparison showing only GNSS and Fused tracks."""
    plt.figure(figsize=(10, 8), dpi=120)
    plt.plot(
        gnss_positions[:, 0],
        gnss_positions[:, 1],
        "-o",
        markersize=3,
        label="GNSS",
        color="blue",
    )
    plt.plot(
        fused_positions[:, 0],
        fused_positions[:, 1],
        "-o",
        markersize=3,
        label="Fused",
        color="green",
    )
    plt.xlabel("East [m]")
    plt.ylabel("North [m]")
    plt.title(title)
    plt.legend()
    plt.axis("equal")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="GNSS + ICP + Fusion Trajectory Tool")
    parser.add_argument(
        "--method",
        choices=["gnss", "icp", "fused", "compare"],
        required=True,
        help="Select trajectory visualization mode",
    )
    parser.add_argument(
        "--max_scans",
        type=int,
        default=None,
        help="Max number of scans to use (default = all)",
    )
    parser.add_argument(
        "--min_distance",
        type=float,
        default=1.0,
        help="Minimum distance (m) to remove LiDAR points near the path in 3D visualization",
    )
    args = parser.parse_args()

    # Output directory
    output_folder = os.path.join("data", "output")
    os.makedirs(output_folder, exist_ok=True)

    # Load PCDs
    pcd_folder = get_default_pcd_folder()
    pcd_dict = load_all_pcds(pcd_folder)

    total_scans = len(pcd_dict)
    max_scans_to_use = args.max_scans or total_scans
    filenames = sorted(pcd_dict.keys())[:max_scans_to_use]

    print(f"Loaded {total_scans} scans. Using {max_scans_to_use}.")

    # GNSS mode
    if args.method == "gnss":
        print("\n=== GNSS-ONLY TRAJECTORY ===")
        odom = GNSSOdometry()
        poses, positions, _ = odom.compute_poses()
        positions = np.array(positions)[:max_scans_to_use]

        # 2D view
        visualize_trajectory_2d(positions, title="GNSS Trajectory (Top-down)")

        # 3D view (with LiDAR frames transformed by GNSS poses)
        points_list = [
            apply_transformation(pcd_dict[fname], poses[i])
            for i, fname in enumerate(filenames)
        ]
        visualize_trajectory_3d(
            positions,
            point_clouds=points_list,
            line_color=(0, 0, 1),
            min_distance=args.min_distance,
        )

        # save merged GNSS-aligned cloud
        output_file = os.path.join(output_folder, "gnss_merged.xyz")
        merged = merge_pointclouds(pcd_folder, poses[:max_scans_to_use])
        save_xyz(merged, output_file)
        print(f"Saved merged GNSS cloud to: {output_file}")

    # ICP mode
    elif args.method == "icp":
        print("\n=== ICP-ONLY TRAJECTORY ===")
        filenames_icp, poses_icp = compute_icp_trajectory(
            pcd_folder, max_scans=max_scans_to_use
        )
        positions_icp = np.array([pose[:3, 3] for pose in poses_icp])

        # 2D view
        visualize_trajectory_2d(positions_icp, title="ICP Trajectory (Top-down)")

        # 3D view (with LiDAR frames transformed by ICP poses)
        points_list = [
            apply_transformation(pcd_dict[fname], poses_icp[i])
            for i, fname in enumerate(filenames_icp)
        ]
        visualize_trajectory_3d(
            positions_icp,
            point_clouds=points_list,
            line_color=(1, 0, 0),
            min_distance=args.min_distance,
        )

        # save merged ICP cloud
        output_file = os.path.join(output_folder, "icp_merged.xyz")
        merged = merge_pointclouds(pcd_folder, poses_icp)
        save_xyz(merged, output_file)
        print(f"Saved merged ICP cloud to: {output_file}")

    # FUSED mode
    elif args.method == "fused":
        print("\n=== FUSED GNSS + ICP TRAJECTORY ===")
        fusion = GNSSICPFusion()
        fused_poses = fusion.fuse(pcd_folder, max_scans=max_scans_to_use)
        fused_positions = np.array([pose[:3, 3] for pose in fused_poses])

        # 2D view
        visualize_trajectory_2d(fused_positions, title="Fused Trajectory (Top-down)")

        # 3D view
        points_list = [
            apply_transformation(pcd_dict[fname], fused_poses[i])
            for i, fname in enumerate(filenames)
        ]
        visualize_trajectory_3d(
            fused_positions,
            point_clouds=points_list,
            line_color=(0, 1, 0),
            min_distance=args.min_distance,
        )

        # save merged fused cloud
        output_file = os.path.join(output_folder, "fused_merged.xyz")
        merged = merge_pointclouds(pcd_folder, fused_poses)
        save_xyz(merged, output_file)
        print(f"Saved merged FUSED cloud to: {output_file}")

    # COMPARE mode: GNSS vs FUSED only
    elif args.method == "compare":
        print("\n=== COMPARING GNSS vs FUSED ===")

        # GNSS
        odom = GNSSOdometry()
        gnss_poses, gnss_positions, _ = odom.compute_poses()
        gnss_positions = np.array(gnss_positions)[:max_scans_to_use]

        # FUSED
        fusion = GNSSICPFusion()
        fused_poses = fusion.fuse(pcd_folder, max_scans=max_scans_to_use)
        fused_positions = np.array([pose[:3, 3] for pose in fused_poses])[
            :max_scans_to_use
        ]

        # Ensure same length
        N = min(len(gnss_positions), len(fused_positions))
        gnss_positions = gnss_positions[:N]
        fused_positions = fused_positions[:N]
        fused_poses = fused_poses[:N]
        filenames_N = filenames[:N]

        # 2D GNSS vs Fused comparison
        visualize_trajectories_2d_comparison(
            gnss_positions, fused_positions, title="GNSS vs Fused (Top-down)"
        )

        # 3D GNSS vs Fused comparison (include fused-transformed LiDAR frames)
        fused_clouds = [
            apply_transformation(pcd_dict[fname], fused_poses[i])
            for i, fname in enumerate(filenames_N)
        ]
        visualize_trajectories_3d_comparison(
            gnss_positions,
            fused_positions,
            point_clouds=fused_clouds,
            min_distance=args.min_distance,
            title="GNSS vs Fused (3D)",
        )

        # Save merged fused cloud (for compare)
        output_file = os.path.join(output_folder, "compare_fused_merged.xyz")
        merged = merge_pointclouds(pcd_folder, fused_poses)
        save_xyz(merged, output_file)
        print(f"Saved merged compare (fused) cloud to: {output_file}")


if __name__ == "__main__":
    main()


# python3 -m src.main --method gnss
# python3 -m src.main --method icp
# python3 -m src.main --method fused
# python3 -m src.main --method compare
