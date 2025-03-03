..
  SPDX-FileCopyrightText: 2019-2024 SeisSol Group

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

* elastic-acoustic simulations are not checkpointed completely (the face displacements are currently ignored)
* the order of the new simulation needs to be identical (this may be subject to change)
