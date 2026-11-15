<!--
SPDX-FileCopyrightText: 2026 SeisSol Group

SPDX-License-Identifier: BSD-3-Clause
SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
-->

# The scripts of the CI

The steps of the jobs in `.github/workflows` are scripts here, so that they
run outside of GitHub as well, e.g. to check a change before pushing it. The
workflows keep what is GitHub's: runners, matrices, checkouts, caches and
artifacts; they hand the values of the matrix to the scripts in the
environment. Each script describes its steps and its environment at its top:

- `jobs/build-variants.sh <arguments of compile_variants.py>`:
  `seissol-build-test` of build-seissol-cpu and build-seissol-gpu
- `jobs/verify.sh`: `seissol-run-test` of both
  (`MPI`, `RANKS`, `THREADS`, `COMMTHREAD`)
- `jobs/perf.sh measure|compare`: `seissol-perf` of build-seissol-cpu
  (`REPETITIONS`)
- `jobs/sanitizer.sh build|test`: sanitizers (`CC`, `CXX`, `SANITIZERS`,
  `CXX_FLAGS`, `CONFIGS`, `DR_QUAD_RULE`, `TESTING_GENERATED`, `THREADS`)
- `jobs/tidy.sh cpu|gpu <part>/4`: clang-tidy (`CLANG_TIDY`)
- `jobs/coverage.sh build|coverage`: coverage (`CC`, `CXX`)
- `jobs/aarch64.sh build|test`: `aarch64` of architectures
  (`CPU`, `HOST_ARCH`)
- `jobs/emulated.sh dependencies|easi|build|test`: `emulated` of
  architectures (`TRIPLET`, `QEMU_CPU`, `CC`, `CXX`, `FLAGS`, `HOST_ARCH`,
  `GEMM_TOOLS`, `CONFIGS`)
- `jobs/pytest.sh [<arguments of pytest>]`: pytest-python (`PYTHON`)
- `jobs/macos.sh build|test`: macos
- `jobs/windows.sh hdf5|easi|build|test`: windows (`DEPS`, `PYTHON`,
  `MPI_HINTS`, `MSMPI_BIN`, `HDF5_VERSION`)

`common.sh` holds what they share (the environment of MPI on the runners),
`build-easi.sh` builds easi as the jobs need it. The checks of the code itself
(clang-format, the headers, flake8, ...) run with
[pre-commit](https://pre-commit.com/): `pre-commit run --all-files`.

## Running a job locally

A script runs where the dependencies of its job are, from any directory; it
works in the checkout it belongs to. E.g. the job of the address and
undefined behavior sanitizers, as its matrix sets it up:

```bash
export CC=gcc CXX=g++ SANITIZERS="address;undefined" \
  CXX_FLAGS=-D_GLIBCXX_ASSERTIONS CONFIGS=ci-sanitize \
  DR_QUAD_RULE=dunavant TESTING_GENERATED=ON THREADS=false
scripts/ci/jobs/sanitizer.sh build && scripts/ci/jobs/sanitizer.sh test
```

The builds of the CPU and the GPU run in the images of the CI
(`seissol/gha-cpu`, `seissol/gha-gpu-nv`, ...; see `images/` here and the
workflows), which have their dependencies, e.g.:

```bash
docker run --rm -v "$PWD:/src" -w /src \
  seissol/gha-cpu:davschneller-gpu-image \
  scripts/ci/jobs/build-variants.sh --preset ci-cpu --toolchain gcc-13 \
  --report build-report.json
```

The other jobs install their dependencies on the runner (see the steps before
the scripts in their workflows); on Ubuntu, the same packages do.

[act](https://github.com/nektos/act) runs the workflows themselves in
Docker, e.g. `act push -W .github/workflows/sanitizers.yml -j sanitizer`
`--matrix setup.name:thread`. The jobs with a container take its image; the
others need one for their runner (`-P ubuntu-26.04=<image>`). The jobs of
macOS and Windows, the ones on self-hosted runners and the ones on a GPU do
not run there. (act has not been tried with these workflows yet.)

## Why a job failed

Failed jobs upload what tells why: the build reports, the results of the
tests (JUnit) with their output, the logs of the builds, the findings of
clang-tidy. GitHub shows the logs of the jobs only after logging in, but
their artifacts without, which `status.py` reads:

```bash
scripts/ci/status.py davschneller/simcopy  # or a commit, or a pull request
```

It prints the checks of the commit, and for each failed job the step that
failed, its failed tests with the first lines of their output that tell what
went wrong, and the findings of clang-tidy. Without a token, GitHub answers
60 requests an hour (a call takes one or two per workflow run with a failed
job); with a token in `GITHUB_TOKEN`, it answers more, and `--logs` prints the
end of the log of each failed job.
