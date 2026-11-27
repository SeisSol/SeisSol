#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The measurement of the kernels with the proxy of an installed SeisSol (seissol-perf in
# .github/workflows/build-seissol-cpu.yml), with scripts/perf/benchmark.py, in the current
# directory:
#
#   scripts/ci/jobs/perf.sh measure   measures into perf/current.json
#   scripts/ci/jobs/perf.sh compare   compares it with perf/baseline.json, if there is one
#
#   INSTALL_DIR   the installation of SeisSol (default: install)
#   REPETITIONS   how often to measure each kernel (default: 5)
#   SUMMARY       where the comparison goes (default: the summary of the step on GitHub, else the
#                 standard output)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"

install_dir="$(realpath "${INSTALL_DIR:-install}")"
benchmark="${SEISSOL_CI_ROOT}/scripts/perf/benchmark.py"
summary="${SUMMARY:-${GITHUB_STEP_SUMMARY:-/dev/stdout}}"

case "${1:-}" in
measure)
  chmod +x "${install_dir}"/bin/*
  export LD_LIBRARY_PATH="${install_dir}/lib:${install_dir}/lib64:${LD_LIBRARY_PATH:-}"
  mkdir -p perf
  python3 "${benchmark}" run \
    --proxy "${install_dir}/bin/proxyseissol" \
    --cells 20000 \
    --timesteps 10 \
    --repetitions "${REPETITIONS:-5}" \
    --warmup 1 \
    --output perf/current.json \
    --github-benchmark perf/github-benchmark.json
  python3 "${benchmark}" noise perf/current.json
  ;;
compare)
  if [ -f perf/baseline.json ]; then
    python3 "${benchmark}" compare \
      --baseline perf/baseline.json \
      --current perf/current.json \
      --summary "${summary}"
  else
    echo "No measurement of master on ${RUNNER_NAME:-this machine} to compare with." >> "${summary}"
  fi
  ;;
*)
  echo "usage: $0 measure|compare" >&2
  exit 2
  ;;
esac
