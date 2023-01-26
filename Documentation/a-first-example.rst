.. _a_first_example:

A first example
===============

This tutorial will guide you through the steps of your first SeisSol
simulation. We will use the `SCEC TPV33
benchmark <http://scecdata.usc.edu/cvws/tpv33docs.html>`__ as an example
in this tutorial. We assume that you have successfully compiled SeisSol.

Setup
-----

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
``OMP_NUM_THREADS=<threads> mpiexec  -np <n> ./SeisSol_<configuration> parameters.par``,
where:

-  ``<configuration>`` depends on your compilation setting (e.g.
   SeisSol_Release_dhsw_4_elastic for a Haswell architecture and order 4 accuracy in space and time).
-  ``<n>`` is the number of MPI ranks / the number of compute nodes used.
-  ``<threads>`` is the number of OpenMP threads per MPI rank, typically the number of CPUS.
   (If you compiled SeisSol with :code:`-DCOMMTHREAD=ON` use the number of CPUs - 1, to reserve one CPU for communication).

**Hint:** Depending on the system you are using, the MPI launcher might
be different from ``mpiexec`` (e.g. ``mpiexec.hydra``, ``mpirun``, ``srun``).
For more infos about how to get optimal performance, have a look at the :ref:`optimal_environment_variables_on_supermuc_ng`.

Result verification
-------------------

SeisSol produces various output files:

* :ref:`3D wave field output <wave_field_output>` (:code:`.xdmf`)
* :ref:`2D free surface output <free_surface_output>` (:code:`-surface.xdmf`)
* :ref:`2D fault output <paraview_output>` (:code:`-fault.xdmf`)
* :ref:`off_fault_receivers` (:code:`-receiver-<id>.dat`)
* :ref:`Fault receivers <fault_receivers>` (:code:`-faultreceiver-<id>.dat`)
* :ref:`energy_output` (:code:`-energy.csv`)

The :code:`xdmf` files can be visualized with `Paraview <https://www.paraview.org/>`__.
For the :code:`dat` files, you can use `viewrec <https://github.com/SeisSol/SeisSol/blob/master/postprocessing/visualization/receiver/bin/viewrec>`__.

The outputs of your simulation can be compared with our outputs (using SeisSol) and the outputs of other codes by checking out the uploaded files for this SCEC benchmark on the SCEC Code Verification Project `website <http://scecdata.usc.edu/cvws/cgi-bin/cvws.cgi>`__.
