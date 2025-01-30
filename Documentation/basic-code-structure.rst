..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Basic code structure
====================

src/
----

============= =============
Folder        Description
============= =============
Common        Routines which are used in several parts of the code, or are generic enough to not fit anywhere else.
Equations     Model-specific code.
Geometry      Everything related to reading tetrahedral meshes and setting up geometry information.
Initializer   Code that is called during initialization, e.g. allocating memory, setting up matrices, parsing material information.
Kernels       Common kernel code.
Model         Common model code.
Modules       Modules system implementation which allows adding code at pre-defined hooks.
Monitoring    Contains code for HPC statistics collected during a run.
Numerical     Helper code for numerics, e.g. quadrature rules.
Parallel      MPI communicator related code.
Physics       Contains friction laws.
Reader        Code for reading parameter files.
ResultWriter  Fault, element, and surface output.
Solver        Time-stepping and code executed during a simulation.
SourceTerm    Everything related to kinematic rupture models.
tests         Unit tests.
============= =============

auto_tuning/
------------

============= =============
Folder        Description
============= =============
config        Sparsity configurations for SeisSol.
proxy         The code directory for SeisSol Proxy. Compiled when compiling SeisSol.
scripts       Wrapper scripts for the SeisSol Proxy.
============= =============

preprocessing/
--------------

============= =============
Folder        Description
============= =============
meshing       Cube generator; gmsh converter; various scripts.
science       ASAGI converter; standard rupture format converter; various scripts.
============= =============

postprocessing/
---------------

============= =============
Folder        Description
============= =============
science       Various scripts processing simulation output.
validation    Cube mesh validation.
visualisation Receiver viewer; scripts.
============= =============
