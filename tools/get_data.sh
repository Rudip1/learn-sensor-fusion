#!/usr/bin/env bash
# Larger datasets are not stored in this repository. This script says where to get them and checks that they
# are in place. Usage:  tools/get_data.sh [nuscenes-mini]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TARGET="${SENSOR_FUSION_EXTERNAL:-$ROOT/data/external}"

case "${1:-nuscenes-mini}" in
  nuscenes-mini)
    DIR="$TARGET/nuscenes"
    if [ -f "$DIR/v1.0-mini/scene.json" ]; then
      echo "nuScenes-mini found in $DIR"
      exit 0
    fi
    cat <<MSG
nuScenes-mini is not in $DIR.

nuScenes may only be downloaded by registered users who accept its terms of use (non-commercial licence):
  1. register and log in at https://www.nuscenes.org/nuscenes#download
  2. download "Mini" (v1.0-mini.tgz) of the full dataset
  3. extract it:   mkdir -p "$DIR" && tar -xzf v1.0-mini.tgz -C "$DIR"
  4. run this script again to check.

The data/external/ folder is ignored by git; never commit the dataset.
MSG
    exit 1
    ;;
  *)
    echo "unknown dataset '$1' (known: nuscenes-mini)" >&2
    exit 2
    ;;
esac
