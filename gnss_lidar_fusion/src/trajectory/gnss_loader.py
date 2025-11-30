# src/trajectory/gnss_loader.py
from typing import Optional
import os
import pandas as pd


class GNSSLoader:
    """
    GNSS CSV loader. Accepts either an explicit path or uses the default dataset path.
    """

    @staticmethod
    def default_csv_path() -> str:
        return os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), "../../data/20251017/parkolo1/fix.csv"
            )
        )

    def load(self, csv_path: Optional[str] = None) -> pd.DataFrame:
        """
        Load CSV and return pandas DataFrame.
        """
        csv_path = csv_path or self.default_csv_path()
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"GNSS CSV not found: {csv_path}")
        # auto-detect separator (tab or comma)
        df = pd.read_csv(csv_path, sep=None, engine="python")
        return df
