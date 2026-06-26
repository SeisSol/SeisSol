..
  SPDX-FileCopyrightText: 2019 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _Checkpointing:

Checkpointing (Beta)
====================

The checkpointing mechanism enables the storage and restart of a simulation to a file. There are (at least) two reasons on why that can be useful:

* fault tolerance. To keep the simulation running after the simulation has failed (or e.g. a node stopped working).
* restart with different parameters

The checkpointing mechanism will automatically adapt to any given element partition. In particular, it is independent of the number of nodes.

Note that, currently, we cannot give any guarantee that the checkpoints are portable between SeisSol versions; use at your own descretion.

Checkpoint Writing
~~~~~~~~~~~~~~~~~~

The checkpoint writing is configured as any other output writer, i.e. you will have to add the ``output`` section the following parameter(s):

.. code:: fortran

   checkpoint = 1
   checkpointinterval = 0.4

Thus, the checkpoint writer is also influenced by the ASYNC parameters.

The checkpoint writer will, in essence, dump the current simulation data into an Hdf5 file, without altering it.

Checkpoint Restart
~~~~~~~~~~~~~~~~~~

The restart is done by adding the command line parameter ``-c <FILE>``. Note that you do not need to activate the checkpointing in the parameter file (i.e. the checkpoint writing as shown above).

The checkpoint reading is then done automatically when starting SeisSol.
Note that the checkpoint reading may take a bit of time: the data needs to be re-distributed, since the element partition will have changed compared to the previous simulation.
Also, the checkpoint will override any set initial condition.

Current Quirks and Limitations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The new checkpointing is still in its beta, hence some interesting behavior is currently allowed.

Generally, a restart is possible if all necessary data is contained in the checkpoint. In particular, it is also possible to ignore data that has been written into a checkpointing file.
The following restart configurations (besides restarting with the same configuration) are therefore possible:

* isotropic/anisotropic elastic to isotropic/anisotropic elastic
* viscoelastic to isotropic/anisotropic elastic (_not_ the other way round)

For the friction laws, the following is the case:

* 6 to 16,1058,33,34
* 16,1058,33,34 to 16,1058,33,34
* 3,4,103 to 3,4,103
* 3,4,103 with Thermal Pressurization to 3,4,103 without Thermal Pressurization


Note however, that the following restrictions apply.

* the order of the new simulation needs to be identical (this may be subject to change)
* the array padding needs to match. (i.e. the checkpoints are architecture-dependent; transfer from CPU to GPU is not supported in general)

Manually Defining a Checkpoint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Last date of update: 2026-06-26.

Note: we do not recommend this procedure, as the format is still subject to change and might be incompatible with future versions of SeisSol. Use with care.

General format outline:

* the checkpoint is a Hdf5 file, stored in the ``checkpoint`` group.
* the group contains two attributes:

  * ``__order``, which is the order that checkpoint was recorded in
  * ``__time``, denoting the simulation time that the checkpoint was recorded at

* for each storage data structure within SeisSol, we have one subgroup. Each subgroup contains an ``__ids`` field which maps the data to the respective cell/face. It is computed as follows:

  * ``lts``: the ID equals the global cell ID.
  * ``dynrup``: the ID is ``globalCellID * 4 + faceSide``. (here, the cell and face side are from the *positive* direction)
  * ``surface``: the ID is ``globalCellID * 4 + faceSide``. (here, we explicitly want a cell face — in fact we may need to store data for both sides of a face here)

* All data fields are "dumped" into the file, meaning that no post processing w.r.t. padding or similar is performed. As general guidelines:

  * the dofs in ``lts`` are given as a two-dimensional array. The fastest-running index runs over the basis functions plus padding; the other index runs over the quantities in their usual SeisSol order.
  * almost all ``dynrup`` fields are defined per quadrature points; meaning that you have a one-dimensional data array per cell which contains all data; plus some padding values. The exact number of data fields varies between the friction laws. Look into https://github.com/SeisSol/SeisSol/blob/master/src/Memory/Descriptor/DynamicRupture.h and for the ``registerCheckpointVariables`` functions for more details.
  * ``surface`` contains the nodal 2D displacements (i.e. again quadrature points) plus padding; they should be mostly relevant for output and/or elastic-acoustic simulations

* The padding is specific to your machine and the precision used; it is usually a power of two (e.g. 8,16,32) and uniform over all data fields. The data sizes are as follows (which are then rounded up by the padding):

  * dofs (order 2 to 8): 4,10,20,35,56,80,120
  * quad points (order 2 to 8; stroud): 9,16,25,36,49,64,81
