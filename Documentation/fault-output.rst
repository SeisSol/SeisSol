.. _fault_output:

Fault output - Dynamic Rupture Visualization Options
============

Introduction
------------

There are two primary methods for visualizing on-fault rupture dynamics:
1. Evaluation of rupture characteristics at specific on-fault locations via ASCII receiver files
2. Visualization across the entire fault surface with files that can be opened with ParaView

While threads or nodes can be allocated to write the ParaView output  (see :ref:`asynchronous-output`), it is typically not required. 

DynamicRupture Namelist OutputPointType Configuration
------------

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
---------------

You can adjust the ParaView output using the Elementwise namelist. Here's an example of how to do this:

.. code-block:: Fortran

  &Elementwise
  printIntervalCriterion = 2 ! 1=iteration, 2=time
  printtimeinterval_sec = 1.0
  OutputMask = 1 1 1 1 1 1 1 1 1 1 1 1 !described herafter
  refinement_strategy = 1 ! or 2
  refinement = 1
  /

printIntervalCriterion
~~~~~~~~~~~~~~~~~~~~~~


- When **printIntervalCriterion** = 1, the output is generated every N time steps, where N is defined by **printInterval**. Note that this option is only compatible with Global time stepping. 
- When **printIntervalCriterion** = 2, output is generated every **printtimeinterval_sec**.

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
11. **DS**: Only with LSW, time at which ASl>D_c (useful for measuring the process zone size)
12. **P_f** and **Tmp**: Only with thermal pressurisation, pore pressure and temperature

SeisSolXdmf Python Module
---------------------

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

.. _ioutputmask-1:

iOutputMask
~~~~~~~~~~~

This is the same as for the ParaView output.

Additional Ascii output
-----------------------

You can output the rupture front at every Gauss point by enabling **RF_output_on** in the DynamicRupture namelist:

.. code-block:: Fortran

  &DynamicRupture
  RF_output_on = 1
  /

We strongly recommend using the ParaView fault output for visualizing the rupture time, as opposed to this ASCII output.
