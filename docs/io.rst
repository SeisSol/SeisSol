..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

The IO Infrastructure (Beta)
============================

The SeisSol IO infrastructure is currently under refactoring, to offer a simpler programming interface, and new IO capabilities.

The new IO infrastructure already covers:

* high-order 3D and 2D output (wavefield, surface, elementwise fault)
* checkpointing of arbitrary friction laws and equations

Support for aggregated point output is under development.

The old IO infrastructure still covers:

* the point output
* zero-th order IO for 3D and 2D output (via the XdmfWriter module)
* metadata like e.g. clustering information per rank

ASYNC
~~~~~

ASYNC serves as the backbone of the new IO infrastructure, and enables the offloading of IO operations to dedicated threads or even processes.
Internally, the new IO system serializes all write requests (to binary or Hdf5 files) and sends them to the ASYNC executors.
See :ref:`asynchronous-output` for more information.
