..
  SPDX-FileCopyrightText: 2018 SeisSol Group

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
IO            Contains the new IO module code. (used for high-order output and checkpointing as of now)
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
============= =============

codegen/
------------

============= =============
Folder        Description
============= =============
config        Sparsity configurations for SeisSol.
kernels       Kernel definitions for Yateto-generated kernels in SeisSol.
matrices      Hard-coded matrix definitions for the Yateto-generated kernels in SeisSol.
============= =============

app/
------------

============= =============
Folder        Description
============= =============
Proxy         Code for the SeisSol Proxy executable. Compiled when compiling SeisSol.
Main          Code for the main SeisSol executable.
============= =============

tests/
------------

Contains unit tests; roughly mirroring the src/ directory structure.

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
