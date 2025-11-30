# src/utils/transforms.py
import numpy as np


def create_transformation_matrix(R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Create 4x4 homogeneous transformation from 3x3 R and 3x1 t.
    """
    T = np.eye(4, dtype=float)
    T[:3, :3] = R
    T[:3, 3] = t.flatten()
    return T


def apply_transformation(points: np.ndarray, T: np.ndarray) -> np.ndarray:
    """
    Apply 4x4 transform to Nx3 points (returns Nx3).
    """
    if points is None or len(points) == 0:
        return points
    n = points.shape[0]
    hom = np.hstack((points, np.ones((n, 1))))
    transformed = (T @ hom.T).T
    return transformed[:, :3]
