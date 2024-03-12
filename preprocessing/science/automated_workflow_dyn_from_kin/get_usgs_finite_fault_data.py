#!/usr/bin/env python3
import json
import os
import wget
import argparse
import shutil


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


def wget_overwrite(url, out_fname=None):
    fn = out_fname if out_fname else os.path.basename(url)
    if os.path.exists(fn):
        os.remove(fn)
    filename = wget.download(url, out=out_fname, bar=None)


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
url = f"https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&&minmagnitude={minM}&eventid={eq_code}"
fn_json = f"{eq_code}.json"
wget_overwrite(url, fn_json)

with open(fn_json) as f:
    jsondata = json.load(f)

properties = find_key_recursive(jsondata, "mag")
for item in properties:
    subjsondata = get_value_by_key(jsondata, item)
    mag = subjsondata
    break

properties = find_key_recursive(jsondata, "place")
for item in properties:
    subjsondata = get_value_by_key(jsondata, item)
    place = subjsondata
    break

properties = find_key_recursive(jsondata, "dyfi")
for item in properties:
    subjsondata = get_value_by_key(jsondata, item)[0]
    eventtime = subjsondata["properties"]["eventtime"]
    break
day = eventtime.split("T")[0]
descr = "_".join(place.split(",")[-1].split())

basic_inversion = find_key_recursive(jsondata, "finite-fault")
for item in basic_inversion:
    subjsondata = get_value_by_key(jsondata, item)[0]
    code_finite_fault = subjsondata["code"]
    update_time = subjsondata["updateTime"]
    hypocenter_x = subjsondata["properties"]["longitude"]
    hypocenter_y = subjsondata["properties"]["latitude"]
    hypocenter_z = subjsondata["properties"]["depth"]

folder_name = f"{day}_Mw{mag}_{descr[:20]}_{code_finite_fault}"

if not os.path.exists(folder_name):
    os.makedirs(folder_name)
if not os.path.exists(f"{folder_name}/tmp"):
    os.makedirs(f"{folder_name}/tmp")


with open(f"{folder_name}/tmp/hypocenter.txt", "w") as f:
    jsondata = f.write(f"{hypocenter_x} {hypocenter_y} {hypocenter_z}\n")

for fn in ["moment_rate.mr", "basic_inversion.param"]:
    url = f"https://earthquake.usgs.gov/product/finite-fault/{code_finite_fault}/us/{update_time}/{fn}"
    wget_overwrite(url, f"{folder_name}/tmp/{fn}")

shutil.move(fn_json, f"{folder_name}/tmp/{fn_json}")
print(folder_name)
