# Data

Small files only (< 5 MB in total). Larger datasets are downloaded by [`tools/get_data.sh`](../tools/get_data.sh)
into ignored folders and never committed.

| path | content | used in |
|---|---|---|
| `gnss/nmea_log1.nmea` | NMEA 0183 log of a consumer GPS receiver (GGA, GSA, GSV, RMC, VTG at 1 Hz), receiver standing still, starting from a cold start, 2025-11-04 | chapters 1, 3 |
| `gnss/nmea_log2.nmea` | same receiver, carried on foot around a loop, 2025-11-04 | chapters 1–3 |
| `parking_lot/run{1,2,3}_fix.csv` | GNSS fixes (ROS `sensor_msgs/NavSatFix`) recorded on a moving platform (TODO: platform type and antenna position) in a parking lot, 2025-10-17: run 1 (165 s, driving), run 2 (8 s, standing), run 3 (106 s, short manoeuvres) | chapters 2, 3, 6 |

The NMEA logs were recorded for this module. They are the receiver's output unchanged except that the
line terminators were normalised to CR LF.

`run*_fix.csv` columns: `stamp_ns` (message time stamp, ns since the Unix epoch), `latitude_deg`,
`longitude_deg`, `altitude_m` (ellipsoidal height, as specified for `NavSatFix`), `status` and `service`
(the `NavSatStatus` fields), and the diagonal of `position_covariance` as `var_east_m2`, `var_north_m2`,
`var_up_m2` (the driver reported a diagonal covariance, `position_covariance_type = 1`). The off-diagonal
entries were zero and the frame ID was `gps` in every message.
