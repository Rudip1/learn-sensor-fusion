import numpy as np
import matplotlib.pyplot as plt
import open3d as o3d
from typing import Optional
from scipy.spatial import cKDTree


# 2D VISUALIZATION (Single Trajectory)
def visualize_trajectory_2d(
    positions: np.ndarray,
    title="Trajectory (Top-down)",
    color="blue",
    linewidth=1,
    linestyle="-",
):
    """Plot a single 2D trajectory (Top-down) with improved styling."""
    plt.figure(figsize=(10, 8), dpi=120)

    plt.plot(
        positions[:, 0],
        positions[:, 1],
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        label=title,
    )

    plt.scatter(
        positions[:, 0],
        positions[:, 1],
        color=color,
        s=2,
        alpha=0.6,
        edgecolor="black",
        linewidth=0.2,
    )

    plt.xlabel("East [m]")
    plt.ylabel("North [m]")
    plt.title(title)
    plt.legend()
    plt.axis("equal")

    plt.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    plt.tight_layout()
    plt.show()


# 2D COMPARISON (GNSS vs FUSED only)
def visualize_trajectories_2d_comparison(
    gnss_positions: np.ndarray,
    fused_positions: np.ndarray,
    title="Trajectory Comparison (GNSS vs Fused)",
):
    plt.figure(figsize=(12, 9), dpi=120)

    def plot_track(pos, color, label):
        plt.plot(pos[:, 0], pos[:, 1], color=color, linewidth=1.5, label=label)
        plt.scatter(pos[:, 0], pos[:, 1], s=3, color=color, alpha=0.3)

    plot_track(gnss_positions, "blue", "GNSS")
    plot_track(fused_positions, "green", "Fused")

    plt.xlabel("East [m]")
    plt.ylabel("North [m]")
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.axis("equal")
    plt.tight_layout()
    plt.show()


# LIDAR POINT REMOVAL (Distance Filter)
def filter_points_near_path(points: np.ndarray, path: np.ndarray, min_distance: float):
    """
    Removes LiDAR points within `min_distance` meters from the trajectory path.
    """
    tree = cKDTree(path[:, :3])
    distances, _ = tree.query(points[:, :3], k=1)

    mask = distances > min_distance  # keep only far points
    return points[mask]


# DISTANCE-COLORED POINT CLOUD (Jet Colormap)
def create_distance_colored_point_cloud(
    points: np.ndarray, path: np.ndarray, min_distance: float = 0.0
) -> o3d.geometry.PointCloud:
    """
    Create point cloud colored by distance to the trajectory path.
    Optionally filter points near the path.
    """
    if min_distance > 0:
        points = filter_points_near_path(points, path, min_distance)

    if len(points) == 0:
        return o3d.geometry.PointCloud()

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(points[:, :3])

    tree = cKDTree(path[:, :3])
    distances, _ = tree.query(points[:, :3], k=1)

    distances_norm = (distances - distances.min()) / (distances.ptp() + 1e-6)
    colors = plt.cm.jet(1 - distances_norm)[:, :3]

    pcd.colors = o3d.utility.Vector3dVector(colors)
    return pcd


# 3D PATH LINES
def create_lineset(positions: np.ndarray, color: tuple):
    """
    Create 3D polyline for trajectory.
    """
    lines = [[i, i + 1] for i in range(len(positions) - 1)]
    line_set = o3d.geometry.LineSet(
        points=o3d.utility.Vector3dVector(positions),
        lines=o3d.utility.Vector2iVector(lines),
    )
    line_set.colors = o3d.utility.Vector3dVector([color for _ in lines])
    return line_set


# 3D SINGLE TRAJECTORY VIEW
def visualize_trajectory_3d(
    positions: np.ndarray,
    point_clouds: Optional[list[np.ndarray]] = None,
    line_color: tuple = (1, 0, 0),
    min_distance: float = 0.0,
):
    """
    Visualize a single 3D trajectory with distance-colored LiDAR.
    """
    geometries = [create_lineset(positions, line_color)]

    if point_clouds:
        for pc in point_clouds:
            geometries.append(
                create_distance_colored_point_cloud(pc, positions, min_distance)
            )

    o3d.visualization.draw_geometries(
        geometries,
        zoom=0.5,
        front=[0.5, 0.5, -1.0],
        lookat=positions[len(positions) // 2],
        up=[0, 0, 1],
    )


# 3D COMPARISON (GNSS vs FUSED only)
def visualize_trajectories_3d_comparison(
    gnss_positions: np.ndarray,
    fused_positions: np.ndarray,
    point_clouds: Optional[list[np.ndarray]] = None,
    title="3D Trajectory Comparison (GNSS vs Fused)",
    min_distance: float = 0.0,
):
    """
    Compare GNSS and Fused trajectories in 3D.
    """
    geometries = []

    # GNSS (blue)
    geometries.append(create_lineset(gnss_positions, (0, 0, 1)))

    # Fused (green)
    geometries.append(create_lineset(fused_positions, (0, 1, 0)))

    # LiDAR
    if point_clouds:
        for pc in point_clouds:
            geometries.append(
                create_distance_colored_point_cloud(pc, fused_positions, min_distance)
            )

    print(f"Opening 3D viewer: {title}")

    o3d.visualization.draw_geometries(
        geometries,
        zoom=0.6,
        front=[0.5, -0.5, -1.0],
        up=[0, 0, 1],
        lookat=fused_positions[len(fused_positions) // 2],
    )
