Basic code structure
====================

src/
----

============= =============
Folder        Description
============= =============
Checkpoint    Code related to checkpointing implementation, which allows to restart a simulation after failure.
Equations     Model-specific code.
Geometry      Everything related to reading tetrahedral meshes and setting up geometry information.
Initializer   Code that is called during initialization, e.g. allocating memory, setting up matrices, parsing material information.
Kernels       Common kernel code.
Model         Common model code.
Modules       Modules system implementation which allows to add code at pre-defined hooks.
Monitoring    Contains code for HPC statistics collected during a run.
Numerical_aux Helper code for numerics, e.g. quadrature rules.
Parallel      MPI communicator related code.
Physics       Contains friction laws.
Reader        Code for reading parameter files.
ResultWriter  Fault, element, and surface output.
Solver        Time-stepping and code executed during a simulation.
SourceTerm    Everything related to kinematic rupture models.
tests         Unit tests.
============= =============

preprocessing/
--------------

============= =============
Folder        Description
============= =============
meshing       Cube generator; gmsh converter; various scripts.
partitioning  *deprecated*
science       ASAGI converter; standard rupture format converter; various scripts.
workflow      *deprecated*
============= =============

postprocessing/
---------------

============= =============
Folder        Description
============= =============
performance   *deprecated*
science       Various scripts processing simulation output.
validation    Cube mesh validation.
visualisation Receiver viewer; scripts.
============= =============
