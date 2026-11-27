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
``fl``          the friction law of the fault faces (kernel mode only)
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

Friction laws
~~~~~~~~~~~~~

The fault faces have the friction law of ``--friction-law``, by its number in the parameter file, with ``tp`` after it for thermal pressurization: one of ``0`` (no friction, the default), ``3``, ``4``, ``6``, ``7``, ``16``, ``33``, ``34``, ``35``, ``103``, ``1058``, and ``3tp``, ``4tp``, ``7tp``, ``103tp``.
All faces have the same parameters, those of a benchmark of the law, under a normal stress of 120 MPa (the parameters of the friction law of the parameter file, and those that the fault files give per point; from these, the initializer of the law derives its initial state, as in a simulation):

=============================== ====================================================================================
Friction law                    Parameters
=============================== ====================================================================================
linear slip weakening (16)      TPV5: ``mu_s`` 0.677, ``mu_d`` 0.525, ``d_c`` 0.4 m, no cohesion; shear stress 81.6 MPa,
                                above the static strength (as in the nucleation patch of TPV5); also for 1058, and
                                with ``pc_vstar`` 0.1, ``pc_prakashLength`` 0.01 for 6
aging and slip law (3, 4)       TPV101: ``rs_a`` 0.008, ``rs_sl0`` 0.02 m, ``rs_f0`` 0.6, ``rs_b`` 0.012, ``rs_sr0``
                                1e-6 m/s, an initial slip rate of 1e-16 m/s at a shear stress of 75 MPa, which the
                                nucleation of TPV101 then raises by 25 MPa on all faces
severe velocity weakening (7)   as the aging law, with ``rs_sl0`` 1e-7 m and a shear stress of 70.8 MPa (+25 MPa)
fast velocity weakening (103)   TPV104: ``rs_a`` 0.01, ``rs_sl0`` 0.4 m, ``rs_srW`` 0.1 m/s, ``rs_muw`` 0.2, ``rs_b``
                                0.014, a shear stress of 40 MPa (+45 MPa)
thermal pressurization (tp)     TPV105: ``tp_halfWidthShearZone`` 0.02 m, ``tp_hydraulicDiffusivity`` 4e-4 m²/s,
                                ``tp_thermalDiffusivity`` 1e-6 m²/s, ``tp_heatCapacity`` 2.7e6, ``tp_undrainedTPResponse``
                                0.1e6, ``tp_iniTemp`` 483.15 K, no initial pore pressure
imposed slip rates (33, 34, 35) 1 m of slip; ``tau_S`` 0.05 s and ``tau_R`` 1 s, or a rise time of 1 s, which began
                                half a rise time before the steps (the delta of 35 comes at the start)
=============================== ====================================================================================

The material on both sides of the faces is that of these benchmarks (2670 kg/m³, 6000 m/s, 3464 m/s), and a time step of 1 ms.
Since the faces read the time derivatives of cells of their own, which the steps do not change, the faults evolve under their friction law alone, without the waves of the cells.
