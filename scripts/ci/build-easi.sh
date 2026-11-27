#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# Builds and installs easi as the CI jobs need it (without ASAGI, Lua and ImpalaJIT), from a checkout
# of it:
#
#   scripts/ci/build-easi.sh <easi source> <prefix> [<further CMake arguments>...]
#
# The build directory is <easi source>/build, unless EASI_BUILD_DIR names another one. The prefix
# has to be writable.

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

source_dir="$1"
prefix="$2"
shift 2
build_dir="${EASI_BUILD_DIR:-${source_dir}/build}"

cmake -S "${source_dir}" -B "${build_dir}" -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${prefix}" \
  -DASAGI=OFF -DLUA=OFF -DEASICUBE=OFF -DIMPALAJIT=OFF \
  "$@"
cmake --build "${build_dir}" --parallel "$(ci_cpus)"
cmake --install "${build_dir}"
