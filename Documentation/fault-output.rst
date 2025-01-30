..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _fault_output:

Fault output
============

Introduction
------------

There are two primary methods for visualizing on-fault rupture dynamics:

1. Evaluation of rupture characteristics at specific on-fault locations via ASCII receiver files
2. Visualization across the entire fault surface with files that can be opened with ParaView

While threads or nodes can be allocated to write the ParaView output  (see :ref:`asynchronous-output`), it is typically not required.

DynamicRupture Namelist OutputPointType Configuration
-----------------------------------------------------

The output type generated is determined by the **OutputPointType** variable in the DynamicRupture namelist. Here is an example of how to configure it:

.. code-block:: Fortran

  &DynamicRupture
  OutputPointType = 4
  /

**OutputPointType** determines the type of output:

| 0 : No output
| 3 : ASCII fault receivers
| 4 : ParaView file
| 5 : Both ASCII fault receivers and ParaView file

.. _paraview_output:

Configuring ParaView output
---------------------------

You can adjust the ParaView output using the Elementwise namelist. Here is an example of how to do this:

.. code-block:: Fortran

  &Elementwise
  printtimeinterval_sec = 1.0
  OutputMask = 1 1 1 1 1 1 1 1 1 1 1 1 !described herafter
  refinement_strategy = 1 ! or 2
  refinement = 1
  vtkorder = -1
  /

printTimeInterval_Sec
~~~~~~~~~~~~~~~~~~~~~

- Output is generated every **printtimeinterval_sec**.

refinement
~~~~~~~~~~

The **refinement** variable determines how the output mesh is created:

- refinement = 0 outputs a single triangle for each mesh cell. The unknowns are calculated at the center of each cell.
- refinement = 1 subdivides each triangle into 3 or 4 subtriangles, depending on the refinement_strategy. A higher refinement value will further subdivide each subtriangle.

- refinement_strategy = 1 divides each triangle into 3 triangles, all sharing the triangle barycenter as a node.
- refinement_strategy = 2 divides each triangle into 4 triangles.

OutputMask
~~~~~~~~~~~

The **OutputMask** variable allows you to select specific unknowns for visualization. You can toggle the output writing of each unknown by changing its corresponding bit in the OutputMask array.

Here's what each bit in the OutputMask array represents:

1. **SRs** and **SRd**: Slip rates in along-strike and along-dip directions
2. **T_s**, **T_d**: Shear stress in strike and dip directions, **P_n**: Normal stress
3. **u_n**: Fault normal velocity (Note: SeisSol does not allow for fault opening (mode I))
4. **Mud**: Current effective friction coefficient, **StV**: State variable for RS friction
5. **Ts0**,\ **Td0**,\ **Pn0**: Total shear and normal stresses, including initial stresses
6. **Sls** and **Sld**: Fault slip in along-strike and -dip directions
7. **Vr**: Rupture velocity, computed from the spatial derivatives of the rupture time
8. **ASl**: Accumulated slip
9. **PSR**: Peak slip rate
10. **RT**: Rupture time
11. **DS**: Dynamic stress time. With LSW, the time at which ASl>D_c. With RS, the time at which mu <= (f0 + mu_w). DS can be used to evaluate the process zone size.
12. **P_f** and **Tmp**: Only with thermal pressurization, pore pressure and temperature

Initial fault tractions
-----------------------

It is worth noticing that **Ts0**,  **Td0**, and  **Pn0** outputted at t=0s are the initial tractions after a first step through the friction routines.
In particular, if the shear traction read in the input files locally exceeds fault strength, the outputted traction at t=0s is reduced compared with the one read in the input files.
To visualize the initial shear tractions (or any other parameters, e.g. d_c) given in the easi file, the script `read_ini_fault_parameter.py <https://github.com/SeisSol/SeisSol/blob/master/preprocessing/science/read_ini_fault_parameter.py>`__ can be used. It depends on the `easi python bindings <https://easyinit.readthedocs.io/en/latest/python_bindings.html>`__.

.. code-block:: bash

    ./read_ini_fault_parameter.py output/data-fault.xdmf fault.yaml --ref_vector " -0.1,0,-1.0"



seissolxdmf python module
-------------------------

You can read SeisSol ParaView files (XDMF/Hdf5 or XDMF/binary files, describing the fault outputs and the free-surface outputs and the volume wavefield outputs) using our Python module **seissolxdmf**. Find it on PyPi at: `seissolxdmf <https://pypi.org/project/seissolxdmf/>`__.

.. _fault_receivers:

Ascii fault receivers
---------------------

To generate ASCII receiver files, configure the **Pickpoint** namelist as in this example:

.. code-block:: Fortran

  &Pickpoint
  printtimeinterval = 1
  OutputMask = 1 1 1 1 1 1 1 1 1 1 1 1 !described herafter
  nOutpoints = 24
  PPFileName = 'fault_receivers.dat'
  /

**printtimeinterval** determines how frequently the output is generated — every **printtimeinterval** (local) time step. Please note that using this output with local time-stepping may result in differently sampled receiver files.

.. _outputmask-1:

OutputMask
~~~~~~~~~~

This is the same as for the ParaView output.

High-Order VTKHDF Output
~~~~~~~~~~~~~~~~~~~~~~~~

The high-order elementwise output can be enabled by setting ``vtkorder`` in the ``elementwise`` section to a positive value, corresponding to the order of the output polynomial per cell.
