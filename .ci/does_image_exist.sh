#!/bin/bash
# SPDX-FileCopyrightText: 2022-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

show_help() {
  echo "Usage - does_image_exist.sh ubuntu 18.04"
}

if [[ -z $1 ]]; then
    echo "ERROR: provide an image name with a tag"
    show_help
    exit 0
fi
image=$1

set +e
docker manifest inspect ${image} > /dev/null

if [ $(echo $?) = '0' ]
then
  echo 1
else
  echo 0
fi
set -e
