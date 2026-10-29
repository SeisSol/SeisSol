<!--
    SPDX-FileCopyrightText: 2021 SeisSol Group

    SPDX-License-Identifier: BSD-3-Clause
    SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

    SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
-->

# Contributing to SeisSol

Thank you for your interest in contributing to SeisSol! Whether you are
fixing a bug, adding a feature, or improving documentation — your help is
appreciated.

This guide explains how to contribute effectively and what to expect during
the review process.

Note that we follow a [Code of Conduct](CODE_OF_CONDUCT.md).
Please be respectful and constructive in all interactions.

## Table of Contents

- [Reporting Bugs](#reporting-bugs)
- [Suggesting Features](#suggesting-features)
- [Contributing to the Code](#contributing-to-the-code)
- [Getting Help](#getting-help)

## Reporting Bugs

When filing a bug report, please include:

- **SeisSol version** (release tag or commit hash — check the first lines
  of SeisSol's output)
- **Build configuration** (compiler, the `CMakeCache.txt` file or
  CMake flags such as `ORDER`, `EQUATIONS`, `PRECISION`,
  `DEVICE_BACKEND`, `HOST_ARCH`)
- **System environment** (OS, HPC cluster name, loaded modules, MPI
  implementation)
- **Steps to reproduce** the issue (parameter file, mesh, easi
  configuration — or a minimal example)
- **Observed vs. expected behavior**
- **Console output or log excerpts** showing errors or warnings

If you suspect a numerical issue, please mention whether the problem
occurs in single precision, double precision, or both.

## Suggesting Features

Feature requests are welcome. Please open an issue and describe:

- The **motivation** or use case behind the feature
- A **proposed solution** (if you have one in mind)
- Any **alternatives** you have considered
- Whether you would be willing to implement it yourself

## Contributing to the Code

### Getting Started

1. **Discuss first.** Before starting significant work, open an
   [issue](https://github.com/SeisSol/SeisSol/issues) or start a
   [discussion](https://github.com/SeisSol/SeisSol/discussions) to
   describe your planned change. This helps us coordinate and avoids
   duplicate effort.

2. **Fork and clone** the repository (**highly important: also include the submodules**):

   ```bash
   git clone --recursive https://github.com/<your-username>/SeisSol.git
   cd SeisSol
   ```

   Or if you have already cloned it, **synchronize**:

   ```bash
   git checkout master
   git pull upstream master
   git submodule update --recursive --init
   ```

3. **Set up your development environment.** Follow the
   [Installing Dependencies](https://seissol.readthedocs.io/en/latest/build-dependencies.html)
   guide, then build SeisSol locally. Make sure you can build
   and run the test suite before making changes.

4. **Install pre-commit hooks**:

   ```bash
   pip install pre-commit
   pre-commit install
   ```

   This ensures formatting and linting checks run automatically before
   each commit. (note: pre-commit might also be available by your package manager)

Now you are ready to join the development.

### Branching Strategy

SeisSol uses a **trunk-based development** model:

- **`master`** is the main development branch. It may be unstable between
  releases.
- **Release tags** (`v1.3.2`, `v1.3.1`, etc.) mark stable versions.
  Bugfix releases are tagged directly from `master` or from short-lived
  release branches when necessary.
- **Feature branches** are created from `master` and merged back via pull
  requests. Name your branch descriptively:
  - `fix/receiver-output-length` for bug fixes
  - `feature/poroelastic-dr` for new features
  - `docs/update-supermuc-guide` for documentation changes
  - `refactor/io-module` for refactoring work

### Working on Your Change

```bash
# Create a feature branch from an up-to-date master
git checkout master
git pull --recurse-submodules
git checkout -b feature/my-new-feature

# Make your changes, then build and test
mkdir -p build && cd build
cmake -DTESTING=ON ..
make -j $(nproc)
ctest --output-on-failure
```

### Code Style

SeisSol consists of the C++ core code (in `src/` and `app/`),
as well as Python code (in `codegen/`) for compile-time
code generation.

#### C++

SeisSol follows C++20, with the following extensions:

- `__restrict` keyword (for better control over array accesses)

C++20 is also used for all CUDA/HIP/SYCL code.

We use `clang-format` and `clang-tidy`
to enforce a consistent code style.

The configuration files (including rule exceptions)
are in [`.clang-format`](.clang-format)
and [`.clang-tidy`](.clang-tidy) at the repository
root.

When submitting code, follow

- **Format your code** before committing: this is done via the pre-commit hooks.
  If you need to format the code manually, run

  ```bash
  clang-format -i src/path/to/your/file.cpp
  ```

- Follow the existing patterns in the codebase and use the `seissol::`
  namespace hierarchy.

- **Prefer modern C++ code**; that is C++20 in our case.

- **Static analysis**: Running `clang-tidy`
  locally before submitting is encouraged, as the CI requires it to
  pass before we can merge your code.
  To enable it, add `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON` to CMake.
  Then run (in the main folder):

  ```bash
  .ci/tidy.sh ./ build/ -fix
  ```

#### Python

Python code (code generation scripts) follows PEP 8,
enforced by `flake8` (configured in [`.flake8`](.flake8),
which also contains a list of currently-excepted rules).
Furthermore, we employ formatting via `black` and `isort`,
and security checks via `bandit`.

All of these tools are run automatically via pre-commit,
and the CI expects them to pass as well.

The scripts in pre/postprocessing can be treated more leniently
(no formatting strictily enforced); but they should ideally
also adhere to these tools.

### Commit Messages

We encourage clear, descriptive commit messages. While we do not strictly
enforce Conventional Commits (nor, admittedly, strictly follow them ourselves),
we recommend the following pattern:

```bash
<type>: <short summary in imperative mood>

<optional body explaining the motivation and context>

<optional footer with references to issues, e.g. "Fixes #1234">
```

Common types: `fix`, `feat`, `docs`, `refactor`, `test`, `perf`, `ci`,
`build`.

Examples:

```bash
fix: correct receiver output length for long simulations

The output buffer was allocated based on the number of time steps
rather than the number of output samples, causing truncation for
simulations with adaptive time stepping.

Fixes #875
```

```bash
feat: add no-fault friction law for GPU backend

Implements the no-fault friction law on the device side, enabling
testing simulations on GPUs without falling back to the host.
```

Keep individual commits focused on a single logical change.
It is fine to have multiple commits in a pull request, but please avoid mixing unrelated
changes (e.g., a bug fix and a formatting change) in the same commit.
But if a sequence of commits concerns the same code change
(i.e. no other commits in between), try combining them
to avoid "inflating" the commit history.
Avoid dedicated formatting commits, unless you port code to a new formatting style.
Extra clang-tidy commits are deemed ok, as it takes longer to apply.

### Pull Request Process

1. **Push your branch** to your fork and open a pull request against
   `SeisSol/SeisSol:master`.

2. **Fill in the PR description** with:
   - A summary of what the PR does and why
   - References to related issues (e.g., "Closes #1234")
   - Any notable design decisions or trade-offs
   - Instructions for testing, if the change is non-trivial
   - AI tools used to create the changes

3. **CI checks** will run automatically. Please ensure:
   - The build succeeds on all CI configurations
   - All existing tests pass
   - Code formatting checks pass (pre-commit (including clang-format, flake8), clang-tidy)

4. **Code review**: At least one maintainer will review your PR. Be
   prepared for feedback — we aim to be constructive and collaborative.
   Reviewers may request changes, ask questions, or suggest alternatives.

5. **Address feedback** by pushing additional commits to your branch.
   Avoid force-pushing during review, as it makes it harder to track
   changes.

6. **Merge**: Once approved, a maintainer will merge your PR. We
   typically use merge commits to preserve the branch history.

### What We Look For in Reviews

- Correctness (does the code do what it claims?)
- Test coverage for new functionality or bug fixes (as far as reasonable)
- Consistency with the existing architecture and code style
- Documentation updates where applicable
- No unnecessary changes to unrelated files

### Testing

#### Running Tests

SeisSol uses CTest for its test suite.
To enable and run tests on an existing CMake installation:

```bash
cmake .. -DTESTING=ON
make -j $(nproc)
ctest --output-on-failure
```

#### Writing Tests

- **Unit tests** for new functionality should be
  placed in the same test subfolder as
  the module they test (e.g., under
  `tests/DynamicRupture/` for `src/DynamicRupture/`).
- Test file names should follow the pattern `Test<ModuleName>.cpp` and `Test<TestSuite>.t.h`.
- If your change fixes a bug, please (if feasible) add a regression test that would
  have caught the issue.
- For changes that affect generated kernels, remember to compare with `-DTESTING_GENERATED=ON`
  in your CMake configuration as a "sanity" check.
  If the generated tests (named "yateto") fail,
  there might be an internal problem in the code
  generation that will need to be treated.

#### Validation

For more comprehensive **end-to-end** validation, the
[Precomputed SeisSol](https://github.com/SeisSol/precomputed-seissol) repository
contains reference setups. If your change affects simulation results,
please verify against relevant benchmarks.

### Code Generation

SeisSol generates optimized matrix kernels at build time using Python
scripts in the `codegen/` directory. If your changes affect the code
generation:

- The generator scripts are invoked by CMake via a custom target
  (`seissol-codegen`). You can re-run the generation step by rebuilding
  this target. In some cases, it is necessary to clean your build directory first.
- Make sure to test with multiple equation types and polynomial orders
  if your change touches the generator.
- Python code generation scripts should also pass `flake8` linting.

### Documentation

SeisSol's documentation lives at
[seissol.readthedocs.io](https://seissol.readthedocs.io) and is built
from Sphinx/reStructuredText files in the `docs/` directory.

If your contribution adds new functionality, a new build parameter, or
changes existing behavior, please update the documentation accordingly.
In particular, if you add a new parameter, add it to `docs/parameters.par` .
You can build the docs locally with:

```bash
cd docs
pip install -r requirements.txt
make html
# Open _build/html/index.html in your browser
```

## Getting Help

If you have questions about contributing or need guidance on where to
start:

- Browse the [Discussions](https://github.com/SeisSol/SeisSol/discussions)
  forum
- Look for issues labeled `good first issue` for beginner-friendly tasks
- Reach out to the maintainers by commenting on an issue

---

Thank you for helping to improve SeisSol!
