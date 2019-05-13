.. _a_first_example:

A first example
===============

This tutorial will guide you through the steps of your first SeisSol
simulation. We will use the `SCEC TPV33
benchmark <http://scecdata.usc.edu/cvws/tpv33docs.html>`__ as an example
in this tutorial. We assume that you have successfully compiled SeisSol
with the options ``parallelization=hybrid`` and ``generatedKernels=yes``
(see :ref:`compiling-seissol`).

Setup
-----

-  Follow the steps 1. - 3. from the [[configuration
   documentation|configuration]].
-  Download the `parameter file and additional setup
   files <https://github.com/SeisSol/Examples/tree/master/tpv33>`__ and
   put them in your launch directory, hereafter named
   ``launch_SeisSol``.
-  Download the mesh `binary
   file <https://syncandshare.lrz.de/getlink/fi72mQiszp6vSs7qN8tdZJf9/tpv33_gmsh>`__
   and the associated `xml
   file <https://syncandshare.lrz.de/getlink/fiEi52Xiwwqkf2sNpTrCHjhw/tpv33_gmsh.xdmf>`__
   and store them in ``launch_SeisSol``.
   **Optional:** For performance reasons, we suggest that you store the
   mesh file in a scratch file system (if one is available at your
   cluster) and create symbolic links in your launch directory (e.g.
   ``ln -s <path/to/tpv33_gmsh> launch_SeisSol/tpv33_gmsh;``). You may
   not see a huge difference in this small test case but for larger
   meshes, this is the recommended strategy.
-  Create the output directory: ``mkdir launch_SeisSol/output``. For the
   output files it might also be beneficial to store them in a scratch
   file system. In this case, create the output directory in your
   scratch file system and a symbolic link ``launch_SeisSol/output`` to
   this directory.

Execution
---------

To execute SeisSol, change to the ``launch_SeisSol`` directory and run:
``OMP_NUM_THREADS=<threads> mpiexec  -np <n> ./SeisSol_<configuration> parameters_<branch>.par``,
where:

-  ``<configuration>`` depends on your compilation setting (e.g.
   SeisSol_release_generatedKernels_dsnb_hybrid_none_9_4 for a Sandy
   Bridge architecture and order 4 accuracy in space and time).
-  ``<n>`` is the number of processes/ the number of partition used.
-  ``<threads>`` is the number of OpenMP threads per process (we usually
   use the number of CPU per core).
-  ``<branch>`` here we provide 2 parameters files, for the master and
   the hardcoded_ini branch. We recommend to use the master branch. The
   master branch relies on easi parameter files (\*.yaml files) for
   setting up some of the simulations properties. On the other hand, the
   hardcoded_ini branch call FORTRAN routines hardcoded in SeisSol to
   set the material and fault properties (in src/Physics/ ini_model.f90
   and ini_model_DR.f90).

**Hint:** Depending on the system you are using, the MPI launcher might
be different from ``mpiexec`` (e.g. ``mpiexec.hydra``).

Result verification
-------------------

The outputs of your simulation can be compared with our outputs (using SeisSol) and the outputs of other codes by checking out the uploaded files for this SCEC benchmark on the SCEC Code Verification Project `website <http://scecdata.usc.edu/cvws/cgi-bin/cvws.cgi>`__.
