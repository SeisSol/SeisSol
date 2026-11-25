#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The clang-tidy jobs (.github/workflows/clang-tidy.yml):
#
#   scripts/ci/jobs/tidy.sh cpu|gpu [<part>/<parts>]
#
# configures build/ as the job of the CPU (ci-cpu) or the GPU (ci-gpu with AdaptiveCpp and
# TensorForge) does, generates the code and checks one part of the translation units with
# scripts/ci/tidy.sh (all of them without a part). Its output goes to build/tidy-output.txt as well,
# its findings to build/tidy-findings.txt.
#
#   CLANG_TIDY          the clang-tidy to run (default: clang-tidy-23)
#   CMAKE_PREFIX_PATH   where the dependencies are (cpu; default: /opt/dependencies)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"
cd "${SEISSOL_CI_ROOT}"

kind="${1:-}"
part="${2:-}"
clang_tidy="$(command -v "${CLANG_TIDY:-clang-tidy-23}")"
"${clang_tidy}" --version

mkdir -p build
cd build
case "${kind}" in
cpu)
  # (the checkout of the job has no submodules)
  git -C .. submodule update --init
  export CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH:-/opt/dependencies}"
  cmake .. -DGEMM_TOOLS_LIST=none -DNETCDF=OFF -DASAGI=OFF -DHDF5=ON -DCMAKE_BUILD_TYPE=Debug \
    -DTESTING=ON -DLOG_LEVEL=warning -DLOG_LEVEL_MASTER=info -DHOST_ARCH=hsw -DCONFIGS=ci-cpu \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DGRAPH_PARTITIONING_LIBS=ptscotch
  make -j "$(ci_cpus)" seissol-codegen
  ;;
gpu)
  CC=gcc CXX=g++ cmake .. -GNinja -DCMAKE_BUILD_TYPE=Debug -DTESTING=ON -DLOG_LEVEL=warning \
    -DLOG_LEVEL_MASTER=info -DHOST_ARCH=hsw -DGEMM_TOOLS_LIST=none -DNETCDF=OFF -DASAGI=OFF -DHDF5=ON \
    -DDEVICE_BACKEND=acpp -DDEVICE_ARCH=sm_75 -DDEVICE_CODEGEN=tensorforge -DCONFIGS=ci-gpu \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
  ninja seissol-codegen
  ;;
*)
  echo "usage: $0 cpu|gpu [<part>/<parts>]" >&2
  exit 2
  ;;
esac

status=0
TIDY_SHARD="${part}" ../scripts/ci/tidy.sh ../ ./ -clang-tidy-binary="${clang_tidy}" -quiet \
  -j "$(ci_cpus)" 2>&1 | tee tidy-output.txt || status=$?
# the findings, without the lines of their context
grep -E "(warning|error): " tidy-output.txt | sort -u > tidy-findings.txt || true
exit "${status}"
