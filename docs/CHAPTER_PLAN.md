# Chapter plan — learn-sensor-fusion   (`<pkg>` = `sensor_fusion`)

GNSS, LiDAR and their fusion for vehicle localization and mapping. Existing material: `gnss_nmea_decoding/`,
`gnss_lidar_fusion/`, `nusc_point_cloud_aggregation/` (move what is worth keeping into the standard layout; the
nuScenes dataset is downloaded by `tools/get_data.sh`, never committed).

1. GNSS fundamentals: how a fix is computed, NMEA 0183 sentences, parsing (C++ parser with tests).
2. Coordinate frames: WGS-84, ECEF, ENU, UTM; conversions with round-trip tests.
3. GNSS quality: DOP, fix types, satellite geometry, detecting motion vs standstill.
4. Point clouds: representation, voxel downsampling, normals, k-d tree (written here, Open3D as reference).
5. Registration: ICP variants, NDT overview; LiDAR odometry and its drift.
6. Loosely coupled GNSS + LiDAR fusion: EKF on SE(2)/SE(3) poses, outlier gating.
7. Map building: aggregating scans with GNSS-INS poses, colouring from cameras (projection).
8. Dynamic objects: ground removal, DBSCAN clustering, temporal consistency filtering.
