#!/usr/bin/env python3
# Copyright (C) 2023 Intel Corporation
# SPDX-License-Identifier: BSD-3-Clause

from argparse import ArgumentParser
import json

parser = ArgumentParser()
parser.add_argument('filenames', nargs='+')
args = parser.parse_args()

for filename in args.filenames:
    with open(filename) as f:
        content = json.load(f)
        print(json.dumps(content, indent=4))
