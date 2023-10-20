# Copyright (C) 2023 Intel Corporation
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required(VERSION 3.13)
project(SeisSol-comm-proxy LANGUAGES CXX)

option(SYCL "Use SYCL to allocate buffers" OFF)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

add_subdirectory(src)
