#!/usr/bin/env python3
input_file_path = "yaml_files/FL33_34_fault.yaml"
output_file_path = "yaml_files/FL33_34_fault_inst.yaml"

with open(input_file_path, "r") as input_file:
    lines = input_file.readlines()

modified_lines = []
for line in lines:
    if 'rupture_onset = x["rupture_onset"]' in line:
        modified_lines.append(
            line.replace('rupture_onset = x["rupture_onset"]', "rupture_onset = 0")
        )
    elif 'rupture_rise_time = x["effective_rise_time"]' in line:
        modified_lines.append(
            line.replace(
                'rupture_rise_time = x["effective_rise_time"]',
                "rupture_rise_time = 0.1",
            )
        )
    else:
        modified_lines.append(line)

with open(output_file_path, "w") as output_file:
    output_file.writelines(modified_lines)
print(f"done writing {output_file_path}")
