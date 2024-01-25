#!/usr/bin/env python3
import json
import os
import wget
import argparse


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
        print(key)
        if isinstance(current_data, dict) and key in current_data:
            current_data = current_data[key]
        else:
            # Key not found, return None or raise an exception based on your needs
            return None

    return current_data


def wget_overwrite(url, out_fname=None):
    print(url)
    fn = out_fname if out_fname else os.path.basename(url)
    if os.path.exists(fn):
        os.remove(fn)
    filename = wget.download(url, out=out_fname)


parser = argparse.ArgumentParser(
    description="download usgs finite model data for a specific earthquake"
)
parser.add_argument("eq_code", help="usgs earthquake code")
parser.add_argument(
    "--min_magnitude", nargs=1, help="min magnitude in eq query", default=[7.0]
)
args = parser.parse_args()


eq_code = args.eq_code
minM = args.min_magnitude[0]
if not os.path.exists(eq_code):
    os.makedirs(eq_code)

url = f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&&minmagnitude={minM}&eventid={eq_code}"
fn_json = f"{eq_code}/{eq_code}.json"
wget_overwrite(url, fn_json)

with open(fn_json) as f:
    jsondata = json.load(f)

# basic_inversion = find_key_recursive(jsondata, 'basic_inversion.param')
basic_inversion = find_key_recursive(jsondata, "finite-fault")
for item in basic_inversion:
    subjsondata = get_value_by_key(jsondata, item)[0]
    code = subjsondata["code"]
    update_time = subjsondata["updateTime"]
    hypocenter_x = subjsondata["properties"]["longitude"]
    hypocenter_y = subjsondata["properties"]["latitude"]
    hypocenter_z = subjsondata["properties"]["depth"]
    with open(f"{eq_code}/hypocenter.txt", "w") as f:
        jsondata = f.write(f"{hypocenter_x} {hypocenter_y} {hypocenter_z}\n")
    print(code, update_time)
    fn = "basic_inversion.param"
    url = (
        f"https://earthquake.usgs.gov/product/finite-fault/{code}/us/{update_time}/{fn}"
    )
    wget_overwrite(url, f"{eq_code}/{fn}")
