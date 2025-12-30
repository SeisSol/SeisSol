..
  SPDX-FileCopyrightText: 2018 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Basic code structure
====================

In general, apply the following rules (later rules supersede earlier rules) when naming files and folders:

- use ``kebab-case`` (i.e. all lowercase; use hyphens to separate words)
- for Python files and directories containing them (e.g. ``codegen/kernels``), use ``snake_case`` (again, all lowercase; but underscores to separate words). That's needed to use these files/folders as modules.
- in C++ source code directories (``src``, ``app``, ``tests``), use ``PascalCase`` (each word starts with a capital letter; but no characters to separate words). That is done due to historic reasons (that is, SeisSol was completely written in FORTRAN in the very beginning).
- if build systems/conventions demand other styles, prefer those. (e.g. ``CMakeLists.txt`` or the ``FindXYZ`` files in the ``cmake`` folder etc.)

src/
----

Keep all folders in ``PascalCase`` here.

=============== =============
Folder          Description
=============== =============
Common          Routines which are used in several parts of the code, or are generic enough to not fit anywhere else.
Dynamic Rupture Contains the Dynamic Rupture implementation (i.e. friction laws, input, output).
Equations       Model-specific code; e.g. flux matrix construction.
Geometry        Everything related to reading tetrahedral meshes and setting up geometry information.
Initializer     Code that is called during initialization, e.g. allocating memory, setting up matrices, parsing material information. Also contains the parameter reader.
IO              Contains the new IO module code. (used for high-order output and checkpointing as of now)
Kernels         Kernel code (i.e. mostly wrappers around the code-generated kernels).
Memory          Helper routines for memory allocation, as well as the storage datastructure and their respective setups.
Model           Common model code.
Modules         Modules system implementation which allows running code at pre-defined hooks.
Monitoring      Contains code for HPC statistics collected during a run.
Numerical       Helper code for numerics, e.g. quadrature rules.
Parallel        MPI/Threading/GPU-related code.
Physics         Hard-coded initial conditions and time inversion.
Reader          Code for reading parameter files.
ResultWriter    The "old" fault, element, surface, and the energy and receiver output.
Solver          Time-stepping and code executed during a simulation.
SourceTerm      Everything related to setting up kinematic rupture models.
=============== =============

codegen/
--------

============= =============
Folder        Description
============= =============
config        Sparsity configurations for SeisSol.
kernels       Kernel definitions for Yateto-generated kernels in SeisSol.
matrices      Hard-coded matrix definitions for the Yateto-generated kernels in SeisSol.
============= =============

app/
----

Keep all folders in ``PascalCase`` here.

============= =============
Folder        Description
============= =============
Proxy         Code for the SeisSol Proxy executable. Compiled when compiling SeisSol.
Main          Code for the main SeisSol executable.
============= =============

tests/
------------

Keep all folders in ``PascalCase`` here.

Contains unit tests; roughly mirroring the src/ directory structure.

docs/
-----

Contains the documentation. Keep all files in ``kebab-case``.

============= =============
Folder        Description
============= =============
doxygen       Doxygen configuration.
figures       Figures and graphics used in the documentation and the readme.
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
