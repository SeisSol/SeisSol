#!/bin/bash

show_help() {
  echo Usage: does_image_exist.sh ubuntu 18.04
}

if [[ -z $1 ]]; then
    echo ERROR: provide a file or a directory
    show_help
    exit 0
fi
path=$1


set +e
git diff --exit-code --quiet HEAD HEAD~1 -- ${path}

if [ $(echo $?) = '0' ]
then
  echo 0
else
  echo 1
fi
set -e
