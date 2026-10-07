import os
import numpy as np


def get_default_pcd_folder():
    return os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../../data/20251017/parkolo1/pcd")
    )


def load_all_pcds(pcd_folder):
    """
    Loads all PCDs as Nx3 numpy arrays.
    Expects .pcd files with simple XYZ (no color).
    """
    pcd_dict = {}
    for f in sorted(os.listdir(pcd_folder)):
        if f.endswith(".pcd"):
            path = os.path.join(pcd_folder, f)
            data = []
            with open(path, "r") as file:
                header_ended = False
                for line in file:
                    if header_ended:
                        parts = line.strip().split()
                        if len(parts) >= 3:
                            data.append(
                                [float(parts[0]), float(parts[1]), float(parts[2])]
                            )
                    elif line.startswith("DATA"):
                        header_ended = True
            if data:
                pcd_dict[f] = np.array(data)
    return pcd_dict
