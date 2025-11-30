import os
import numpy as np
import pandas as pd
from .pcd_loader import load_all_pcds, get_default_pcd_folder


def read_gnss_csv(csv_file):
    return pd.read_csv(csv_file)


def filename_to_timestamp(filename):
    return int(filename.replace(".pcd", ""))


def match_pcd_to_gnss(pcd_filenames, gnss_df):
    gnss_ts = gnss_df["header_stamp_secs"] * 1e9 + gnss_df["header_stamp_nsecs"]
    matched = {}
    for f in pcd_filenames:
        idx = (np.abs(gnss_ts - filename_to_timestamp(f))).argmin()
        matched[f] = idx
    return matched
