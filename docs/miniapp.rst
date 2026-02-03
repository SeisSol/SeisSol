..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

SeisSol Proxy (mini app)
========================

SeisSol has a mini app to test the performance of its kernels.

Internally, SeisSol uses the mini app to estimate the runtime (called "Mini SeisSol" in this context).

It is also built as a standalone application called "SeisSol Proxy" alongside the main SeisSol executable.

SeisSol Proxy Usage
~~~~~~~~~~~~~~~~~~~

The SeisSol Proxy runs a kernel (or more, in sequence) over several cells and repeats this process for a fixed number of times.

In other words, the SeisSol proxy runs a single time cluster (i.e. single loops over all cells);
thus equalling a simulation using a global time stepping scheme.

It supports output (``--format`` parameter) on the command line (format ``plain``, and ``plaintflop``), and in a JSON format (``json``). We recommend using the latter (JSON)
to gather performance data.
Normally, we use ``GFLOP(/s)`` as basis for the output values; unless you rescale to ``TFLOP(/s)`` with ``plaintflop`` (in case ``GFLOP`` becomes too large).

The basic usage parameters are as follows:

.. code-block:: bash

  ./proxyseissol-<yourconfig> <cell count> <repetition count> <kernels>

It should be noted that internally, the proxy will run through the kernels one time *before* it starts
with the actual repetitions, to potentially initialize kernels and data. This runthrough is not
counted into the final time. However, it might be good to keep its existence in mind when profiling kernels.

As an example, the following command will run the ``all`` kernel sequence over 100000 for 10 times (excluding the initialization loop over the 100000 cells).

.. code-block:: bash

  ./proxyseissol-<yourconfig> 100000 10 all

For the list of available kernels, we refer to the command line documentation.
Note that ``all`` is equivalent to running ``local,neigh`` in this case.

Kernel Chaining
~~~~~~~~~~~~~~~

The proxy also supports combining and chaining multiple kernels. For that, write them with a comma between them.
For example, ``local,neigh,neigh_dr`` will run ``local``, then ``neigh``, then ``neigh_dr`` in this sequence.
