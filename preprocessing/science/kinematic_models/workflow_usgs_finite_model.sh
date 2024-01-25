#!/bin/bash

event_id=$1
python get_usgs_finite_fault_data.py $event_id
cd $event_id
read lon lat _ < hypocenter.txt
proj="+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=${lon} +lat_0=${lat}"
echo $proj
../generate_FL33_input_files.py basic_inversion.param --spatial_zoom 3 --generate_ts_yaml "${proj}"
../generate_usgs_finite_fault_mesh.py
