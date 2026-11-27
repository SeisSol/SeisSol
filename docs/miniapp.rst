..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _miniapp:

SeisSol Proxy (mini app)
========================

SeisSol has a mini app to measure the performance of its kernels.
It is built as a standalone application, the SeisSol Proxy (``proxyseissol``, or ``SeisSol_proxy_<configuration>`` without ``NEW_BINARY_NAMING``), alongside the main SeisSol executable.
SeisSol runs it internally as well, to estimate how fast its ranks compute ("Mini SeisSol", see ``SEISSOL_MINISEISSOL`` in :ref:`env_vars`).

Usage
~~~~~

The proxy runs one or more kernels over a number of cells, and repeats this a number of times.
It runs a single time cluster, i.e. loops over all cells as a simulation with global time stepping does, on synthetic data: cells with random neighbors and arbitrary values.

.. code-block:: bash

  ./proxyseissol <cells> <repetitions> <kernels>

Before the repetitions, the proxy runs through the kernels once, to initialize kernels and data; this run does not count into the time.
For example, the following command runs the ``all`` kernels over 100000 cells 10 times:

.. code-block:: bash

  ./proxyseissol 100000 10 all

The proxy writes its output as text (``-f plain``, the default) or as JSON (``-f json``); we recommend the latter for gathering performance data.
``-c <configuration>`` chooses one of the configurations built into the executable (see ``--help``).

Kernels and their chaining
~~~~~~~~~~~~~~~~~~~~~~~~~~

The kernels are:

=============== ==================================================================================
Kernel          Description
=============== ==================================================================================
``local``       the time and the local integration of the cells
``neigh``       the contributions of their neighbors
``neigh_dr``    the same with all faces of the cells on a fault
``godunov_dr``  the work of the fault faces (see the modes below)
``all``         ``local,neigh``: a time step without dynamic rupture
``all_dr``      ``local,godunov_dr,neigh_dr``: a time step with dynamic rupture
``ader``        the time integration of the cells, the first part of ``local`` (kernel mode only)
``localwoader`` the local integration of the cells, the second part of ``local`` (kernel mode only)
``plast``       the plastic correction of cells that stay within the yield surface (kernel mode only)
``plast_yield`` the plastic correction of cells that yield in every step (kernel mode only)
=============== ==================================================================================

Kernels can be chained with commas: ``local,neigh,neigh_dr`` runs ``local``, then ``neigh``, then ``neigh_dr`` in each repetition.
``plast`` and ``plast_yield`` need different cells, and cannot be chained with each other.

Modes
~~~~~

The proxy runs the kernels in one of two modes (``--mode``):

- ``kernel`` (the default): each kernel runs over all cells or faces on its own, with nothing of what a simulation does around it. ``godunov_dr`` is the space-time interpolation of the fault faces alone.
- ``cluster``: the loops of the clusters of a simulation run, with all they do around the kernels, step by step and outside of the time stepping. ``local`` is the prediction of a cluster of cells, ``neigh`` its correction, and ``godunov_dr`` the work of a cluster of fault faces, the friction law included (a fault without friction). With ``--plasticity elastic`` or ``--plasticity yield``, the cells have a plasticity, which their correction runs: cells that stay within the yield surface, or cells that yield in every step.

The two modes count the same FLOPs for the same kernels; the friction law counts none.
Mini SeisSol runs in the cluster mode.
