..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

The IO Infrastructure
=====================

SeisSol internally uses a single IO module to handle all output-related issues.
It enables us to send write requests asynchronously to the file system; and not
disrupt the computation.
In particular, the IO system covers:

* high-order 3D and 2D output (wavefield, surface, elementwise fault)
* zero-th order IO for 3D and 2D output (via the XdmfWriter module)
* checkpointing of arbitrary friction laws and equations

Not directly handled by it are:

* the point output (on-fault and off-fault)
* metadata like e.g. clustering information per rank

ASYNC
~~~~~

ASYNC serves as the backbone of the IO infrastructure, and enables the offloading of IO operations to dedicated threads or even processes.
Internally, the IO system serializes all write requests (to binary or Hdf5 files) and sends them to the ASYNC executors.
See :ref:`asynchronous-output` for more information.
