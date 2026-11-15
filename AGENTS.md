<!--
SPDX-FileCopyrightText: 2026 SeisSol Group

SPDX-License-Identifier: BSD-3-Clause
SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
-->

# Checking a change before it goes to the CI

For everyone who changes SeisSol, people and coding agents alike: what the
CI checks, how to check it before pushing, and where a local check differs
from the CI. [CONTRIBUTING.md](CONTRIBUTING.md) has the rules of the code
itself; [scripts/ci/README.md](scripts/ci/README.md) has the scripts of the
jobs of the CI.

## Before pushing

1. `pre-commit run --all-files`, or the installed hook (`pre-commit
   install`): clang-format, the headers and names of the files, flake8,
   black, isort, markdownlint and more. Commit what it rewrites, too.
2. Build the configurations that the change touches with `TESTING=ON`, and
   run their tests: `ctest -L <label>` or `ctest -R <name>`, e.g.
   `ctest -R "verification\.smoke"`. The CI builds the sets of
   `cmake/confignames.cmake` (`ci-cpu`, `ci-gpu`, ...), each into one
   executable.
3. clang-tidy on the translation units that the change touches (see
   below); the CI runs it on all of them, for the CPU (`ci-cpu`) and for the
   GPU (`ci-gpu` with AdaptiveCpp).
4. A change of Python code: `scripts/ci/jobs/pytest.sh`.

## clang-tidy

`scripts/ci/jobs/tidy.sh cpu|gpu` configures and checks as the CI does. To
check only some translation units, configure as the job does (with
`-DCMAKE_EXPORT_COMPILE_COMMANDS=ON`), generate the code
(`ninja seissol-codegen`), and run clang-tidy on them with the header filter
of the CI:

```bash
clang-tidy -p build --quiet -header-filter='^$|app/|src/|tests/' src/X.cpp
```

The translation units that include a changed header come from
`clang-scan-deps -compilation-database build/compile_commands.json -format
make`. A finding in a header fails every part of the job that reads it.
Newer C libraries than the one of the CI let clang-tidy report
`misc-static-assert` for each `assert(x > 0)`; the CI does not.

## Where a local check differs from the CI

- **The kernels.** The images of the CI have LIBXSMM, whose kernels load
  the tensors that are declared aligned with aligned instructions; a build
  with PSpaMM alone does not. An array that is not aligned crashes there
  only (SIGSEGV at address 0).
- **The verification cases.** `tests/verification/run_case.py` runs a case
  in its configuration (`SEISSOL_CONFIGURATION`); a run by hand without it
  takes the first configuration of the executable. `ctest -R <test> -V -N`
  prints the command line of a case.
- **The snapshot references.** The references of the snapshot cases are
  recorded in one environment and compared bitwise there; the CI compares
  them with a tolerance (`SEISSOL_SNAPSHOT_TOLERANCE_DOUBLE=1e-9`), which a
  local build needs as well. `SEISSOL_UPDATE_SNAPSHOTS=1 ctest -R
  snapshot...` records them anew, with the commit of a clean tree.
- **MPI on a GPU.** Several ranks on a GPU need an MPI that can read device
  memory, or `SEISSOL_TRANSFER_MODE=host`. Open MPI 5 with more than one
  accelerator framework stops in MPI_Init: `OMPI_MCA_accelerator=cuda` (or
  `rocm`) chooses one; `null` makes it unable to read device memory.
- **The runners.** The runners of GitHub have 4 cores and 16 GB: builds
  with many configurations can run out of memory there (NVHPC in
  particular), and tests take longer. Processes started side by side with
  `OMP_PLACES=cores` and one thread all end up on the same core.
- **The pinned memory of GPU builds.** Host memory that ends at a page that
  may not be touched turns a read past the end of an array into a crash,
  which other builds do not see.

## After pushing

`scripts/ci/status.py <branch>` prints the checks of the branch and, for each
failed job, its failed step, failed tests with what went wrong, and the
findings of clang-tidy, without logging in to GitHub. A test that fails in
the CI but passes locally once is no proof of a flaky test: run it again in
the image of the job, under load (`ctest -j`), before calling it that.

## Commits

The history uses `type: summary in lowercase` (`fix`, `feat`, `refactor`,
`test`, `docs`, `ci`, `verification`, ...) with a body in prose that tells
the cause and the effect. Commits made with a coding agent end with a line
that says so, e.g. `AI-generated. Model: <model>`.
