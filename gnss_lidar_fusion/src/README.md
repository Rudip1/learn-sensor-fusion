# GNSS + LiDAR Fusion Project

This project demonstrates GNSS, LiDAR, and fusion-based trajectory estimation and visualization.  
It includes scripts for computing trajectories using GNSS, ICP, and GNSS+ICP fusion, along with point cloud visualization.

---

## **Folder Structure**
```bash
gnss_lidar_fusion/
│
├── data/ # Input data (LiDAR scans, GNSS logs)
│ └── output/ # Outputs (merged point clouds, .xyz files)
├── report/ # Reports, analysis, or plots
├── src/ # Source code
│ ├── main.py # Main script to run trajectories
│ ├── trajectory/ # Trajectory computation modules
│ └── utils/ # Utility functions (PCD loading, transforms, visualization)
├── requirements.txt # Python dependencies
├── run.sh
└── README.md
```

---

## **Setup**

1. Clone the repository:
```bash
git clone git@github.com:Rudip1/gnss_lidar_fusion.git
cd gnss_lidar_fusion
```

2. Create a Python virtual environment (recommended):
```bash
python3 -m venv venv
source venv/bin/activate
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## **Run Trajectories**

1. GNSS-only trajectory
```bash
python3 -m src.main --method gnss
```

2. ICP-only trajectory
```bash
python3 -m src.main --method icp
```

3. Fused GNSS + ICP trajectory
```bash
python3 -m src.main --method fused
```

4. Compare GNSS and fused trajectories
```bash
python3 -m src.main --method compare
```


All merged point clouds will be saved in:
```bash
data/output/
```

Notes
 .xyz merged point clouds in data/output/


## **Outputs**

---

### 1 GNSS-only Trajectory

**2D Top-down view:**  
![GNSS 2D](https://github.com/Rudip1/gnss_lidar_fusion/blob/main/report/gnss_2d.png?raw=true)

**3D Trajectory view:**  
![GNSS 3D](https://github.com/Rudip1/gnss_lidar_fusion/blob/main/report/gnss_3d.png?raw=true)

---

### 2 ICP-only Trajectory

**2D Top-down view:**  
![ICP 2D](https://github.com/Rudip1/gnss_lidar_fusion/blob/main/report/icp_2d.png?raw=true)

**3D Trajectory view:**  
![ICP 3D](https://github.com/Rudip1/gnss_lidar_fusion/blob/main/report/icp_3d.png?raw=true)

---

### 3 Fused GNSS + ICP Trajectory

**2D Top-down view:**  
![Fused 2D](https://github.com/Rudip1/gnss_lidar_fusion/blob/main/report/fused_2d.png?raw=true)

**3D Trajectory view:**  
![Fused 3D](https://github.com/Rudip1/gnss_lidar_fusion/blob/main/report/fused_3d.png?raw=true)

---

### 4 Comparison: GNSS vs Fused Trajectories

**2D Top-down comparison:**  
![Comparison 2D](https://github.com/Rudip1/gnss_lidar_fusion/blob/main/report/compare_2d.png?raw=true)

**3D Trajectory comparison:**  
![Comparison 3D](https://github.com/Rudip1/gnss_lidar_fusion/blob/main/report/compare_3d.png?raw=true)
