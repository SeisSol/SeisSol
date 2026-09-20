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
* zero-th order IO for 3D and 2D output
* checkpointing of arbitrary friction laws and equations

Every mesh output writes one self-contained file per output step by default. Two
other groupings are available through ``outputtimeseries``, or per output through
``wavefieldtimeseries``, ``surfacetimeseries`` and ``timeseries`` in the
``elementwise`` section, and both need the output to be written as VTKHDF (that
is, with a ``vtkorder`` set):

* ``incremental`` writes the data that does not change -- the points and the
  connectivity -- once per run, into a file of its own, and lets every snapshot
  reference it. Worth it whenever the mesh is large against the fields written
  on it.
* ``monolith`` writes one file for the whole run, as a VTKHDF time series. The
  step bookkeeping lives inside that file, so no ``.pvd`` collection is written
  alongside it.

Not directly handled by it are:

* the point output (on-fault and off-fault)
* metadata like e.g. clustering information per rank

ASYNC
~~~~~

ASYNC serves as the backbone of the IO infrastructure, and enables the offloading of IO operations to dedicated threads or even processes.
Internally, the IO system serializes all write requests (to binary or Hdf5 files) and sends them to the ASYNC executors.
See :ref:`asynchronous-output` for more information.
