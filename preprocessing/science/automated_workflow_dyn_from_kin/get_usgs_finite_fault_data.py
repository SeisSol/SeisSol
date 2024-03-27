#!/usr/bin/env python3
import json
import os
import wget
import argparse
import shutil
import numpy as np
from obspy import UTCDateTime


def find_key_recursive(data, target_key, current_path=None):
    if current_path is None:
        current_path = []

    occurrences = []

    if isinstance(data, dict):
        for key, value in data.items():
            new_path = current_path + [str(key)]
            if key == target_key:
                occurrences.append("|".join(new_path))
            occurrences.extend(find_key_recursive(value, target_key, new_path))
    elif isinstance(data, list):
        for index, item in enumerate(data):
            new_path = current_path + [str(index)]
            occurrences.extend(find_key_recursive(item, target_key, new_path))

    return occurrences


def get_value_by_key(data, target_key):
    keys = target_key.split("|")
    current_data = data

    for key in keys:
        if isinstance(current_data, dict) and key in current_data:
            current_data = current_data[key]
        else:
            # Key not found, return None or raise an exception based on your needs
            return None

    return current_data


def get_value_from_usgs_data(jsondata, key):
    item = find_key_recursive(jsondata, key)[0]
    return get_value_by_key(jsondata, item)


def wget_overwrite(url, out_fname=None):
    fn = out_fname if out_fname else os.path.basename(url)
    if os.path.exists(fn):
        os.remove(fn)
    filename = wget.download(url, out=out_fname, bar=None)


def retrieve_usgs_id_from_dtgeo_dict(fname, min_mag):
    ev = np.load(fname, allow_pickle=True).item()

    origin_time = UTCDateTime(ev["ot"])
    origin_time_plus_one_day = origin_time + 24 * 3600

    url = f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson"
    url += f"&minmagnitude={min_mag}"
    url += f"&starttime={origin_time}&endtime={origin_time_plus_one_day}"
    url += f"&minlatitude={ev['lat']-0.5}&maxlatitude={ev['lat']+0.5}"
    url += f"&minlongitude={ev['lon']-0.5}&maxlongitude={ev['lon']+0.5}"

    fn_json = "out.json"
    wget_overwrite(url, fn_json)

    with open(fn_json) as f:
        jsondata = json.load(f)
    features = get_value_from_usgs_data(jsondata, "features")
    if not features:

        def convert(obj):
            # Convert NumPy arrays to lists before pretty printing
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj

        pretty_dict = json.dumps(dict(ev), indent=4, sort_keys=True, default=convert)
        print(pretty_dict)
        raise ValueError(f"usgs event_id could not be retrieved from {fname}")
    usgs_id = features[0]["id"]
    return usgs_id


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="download usgs finite model data for a specific earthquake"
    )
    parser.add_argument(
        "usgs_id_or_dtgeo_npy",
        help="usgs earthquake code or event dictionnary (dtgeo workflow)",
    )
    parser.add_argument(
        "--min_magnitude",
        nargs=1,
        help="min magnitude in eq query",
        default=[7.0],
        type=float,
    )
    parser.add_argument("--suffix", nargs=1, help="suffix for folder name")
    args = parser.parse_args()
    minM = args.min_magnitude[0]

    if args.usgs_id_or_dtgeo_npy[-3:] == "npy":
        usgs_id = retrieve_usgs_id_from_dtgeo_dict(args.usgs_id_or_dtgeo_npy, minM)
    else:
        usgs_id = args.usgs_id_or_dtgeo_npy

    url = f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&&minmagnitude={minM}&eventid={usgs_id}"
    fn_json = f"{usgs_id}.json"
    wget_overwrite(url, fn_json)

    with open(fn_json) as f:
        jsondata = json.load(f)

    mag = get_value_from_usgs_data(jsondata, "mag")
    place = get_value_from_usgs_data(jsondata, "place")
    dyfi = get_value_from_usgs_data(jsondata, "dyfi")[0]
    eventtime = dyfi["properties"]["eventtime"]
    day = eventtime.split("T")[0]
    descr = "_".join(place.split(",")[-1].split())
    finite_fault = get_value_from_usgs_data(jsondata, "finite-fault")[0]
    code_finite_fault = finite_fault["code"]
    update_time = finite_fault["updateTime"]
    hypocenter_x = finite_fault["properties"]["longitude"]
    hypocenter_y = finite_fault["properties"]["latitude"]
    hypocenter_z = finite_fault["properties"]["depth"]
    suffix = args.suffix[0] if args.suffix else ""

    folder_name = f"{day}_Mw{mag}_{descr[:20]}_{code_finite_fault}{suffix}"

    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    if not os.path.exists(f"{folder_name}/tmp"):
        os.makedirs(f"{folder_name}/tmp")

    with open(f"{folder_name}/tmp/hypocenter.txt", "w") as f:
        jsondata = f.write(f"{hypocenter_x} {hypocenter_y} {hypocenter_z}\n")

    for fn in ["moment_rate.mr", "basic_inversion.param", "complete_inversion.fsp"]:
        url = f"https://earthquake.usgs.gov/product/finite-fault/{code_finite_fault}/us/{update_time}/{fn}"
        wget_overwrite(url, f"{folder_name}/tmp/{fn}")

    shutil.move(fn_json, f"{folder_name}/tmp/{fn_json}")
    print(folder_name)
