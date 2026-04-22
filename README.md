<!--
    SPDX-FileCopyrightText: 2012 SeisSol Group

    SPDX-License-Identifier: BSD-3-Clause
    SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

    SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
-->

# ![SeisSol](docs/figures/logo-sans-darkred-border.svg)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4672483.svg)](https://doi.org/10.5281/zenodo.4672483)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![GitHub Repo stars](https://img.shields.io/github/stars/SeisSol/SeisSol)](https://github.com/SeisSol/SeisSol/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/SeisSol/SeisSol)](https://github.com/SeisSol/SeisSol/network/members)
[![GitHub issues](https://img.shields.io/github/issues/SeisSol/SeisSol)](https://github.com/SeisSol/SeisSol/issues)
[![GitHub pull requests](https://img.shields.io/github/issues-pr/SeisSol/SeisSol)](https://github.com/SeisSol/SeisSol/pulls)
[![Documentation status](https://app.readthedocs.org/projects/seissol/badge)](https://seissol.readthedocs.io)

SeisSol is a scientific software for the numerical simulation of seismic wave
phenomena and earthquake dynamics. It is based on the discontinuous Galerkin
method combined with ADER time discretization.

> [!NOTE]
> SeisSol is still under heavy development and comes without any guaranteed
> functionality. At the moment we can only provide very limited support for
> general users.

## Quick Start

### Prerequisites

- C++ compiler (GCC ≥ 9, Clang ≥ 18, ICX, or NVHPC)
- CMake ≥ 3.20
- Python ≥ 3.9 with numpy and setuptools
- MPI (standard ≥ 2.2), HDF5 (parallel), Eigen ≥ 3.4,
and [easi](https://github.com/SeisSol/easi) ≥ 1.5
- (recommended) ParMETIS or another supported graph partitioning library
- (recommended) LIBXSMM, PSpaMM, (for GPUs) TensorForge

For a complete list of dependencies and installation
methods (including Spack), see the
[Installing Dependencies](https://seissol.readthedocs.io/en/latest/build-dependencies.html)
guide.

### Building SeisSol

```bash
# Clone the repository with all submodules (important!)
git clone --recursive https://github.com/SeisSol/SeisSol.git
cd SeisSol

# Configure and build (CPU example, 4th order, elastic)
mkdir -p build && cd build
cmake -DORDER=4 -DEQUATIONS=elastic -DPRECISION=double ..
make -j $(nproc)
```

For GPU builds, add the device backend and
architecture flags (here CUDA, Hopper GPU):

```bash
cmake -DORDER=4 -DPRECISION=single \
      -DDEVICE_BACKEND=cuda -DDEVICE_ARCH=sm_90 ..
```

See the [Build Documentation](https://seissol.readthedocs.io/en/latest/build-overview.html)
for detailed instructions and the
[Build Parameters Reference](https://seissol.readthedocs.io/en/latest/build-parameters.html)
for all available configuration options.

The full documentation [seissol.readthedocs.io](https://seissol.readthedocs.io)
also contains setup guides for specific HPC clusters
(e.g. SuperMUC-NG, Leonardo, Frontera, and others).

## Documentation

Visit [seissol.readthedocs.io](https://seissol.readthedocs.io)
for the full documentation, including a detailed overview
over the meshing process and available features and options.

## Community and Support

- **Discussion Forum**: Ask questions, share ideas, and
connect with developers and users in our [GitHub Discussions](https://github.com/SeisSol/SeisSol/discussions).
- **Bug Reports**: If there is a problem with SeisSol,
please open an [issue](https://github.com/SeisSol/SeisSol/issues)
with a description of the problem, your build configuration, and steps to reproduce.
- **Website**: [www.seissol.org](http://www.seissol.org)

## Contributing

We welcome contributions of new features, extensions, and bug fixes.
Before starting, please discuss your plans by opening an issue so we can
coordinate effectively.

For detailed guidelines on code style, branching, testing, and the pull
request workflow, see [CONTRIBUTING.md](CONTRIBUTING.md).

## Citing SeisSol

If you utilize SeisSol in a publication or want to refer to it,
please follow the suggestions on our [How To Cite](https://seissol.org/about/howtocite/)
page.
It also includes examples how to cite specific SeisSol features,
models or previous use cases.

To reference SeisSol as a software package and the specific version you used,
please provide the link [doi.org/10.5281/zenodo.4672483](https://doi.org/10.5281/zenodo.4672483)
which points to Zenodo.

## Collaboration

If you are interested in a close collaboration, please contact [Alice Gabriel](https://www.alicegabriel.com/).

## Code of Conduct

We follow a [Code of Conduct](CODE_OF_CONDUCT.md). Please adhere to it
when participating in our community.

## Licensing

SeisSol is licensed under the [BSD-3-Clause License](LICENSE).
Some files in the `cmake` and `external` folders may carry different
licenses (BSL-1.0, MIT, BSD-2-Clause). See the individual file headers
and the [LICENSES](LICENSES/) directory for details.
