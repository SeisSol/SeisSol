..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Postprocessing and Visualization
================================

To visualize the output of SeisSol, you can use ParaView with parallel
rendering servers to distribute the workload. This chapter documents how
to achieve this on the SuperMUC cluster. For other clusters, the
approach might be similar.

Login to SuperMUC and start the pvservers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Connect to the SuperMUC via

.. code-block:: bash

   ssh -Y hw.supermuc.lrz.de -l username

If you don't have access please see
`https://www.lrz.de/services/compute/supermuc/access_and_login/ <https://www.lrz.de/services/compute/supermuc/access_and_login/>`__
for support.

To start the pvserver processes on the cluster use a jobscript to
specify all the parameters. Here is an example of how this can look like
(save as start_paraview.sh):

.. code-block:: bash

   #!/bin/bash
   ##
   ###  optional: energy policy tags
   ##
   ##  DO NOT USE environment = COPY_ALL
   #@ job_type = MPICH
   #@ class = test
   #@ node = 2
   ####  schedule the job to exactly 1 island
   ########## (more than one will fail)
   #@ island_count=1
   #@ tasks_per_node = 28
   #@ wall_clock_limit = 0:30:00
   #@ job_name = paraview
   #@ network.MPI = sn_all,not_shared,us
   #@ initialdir = $(home)/viz
   #@ output = job.$(schedd_host).$(jobid).out
   #@ error =  job.$(schedd_host).$(jobid).err
   #@ queue
   . /etc/profile
   . /etc/profile.d/modules.sh
   ## setup of environment
   module unload mpi.ibm
   module load mpi.intel paraview/5.2.0_mesa
   mpiexec -n 56 pvserver --use-offscreen-rendering

It is important to use a ParaView version that has been built with MPI
and MESA support because usually the compute nodes in the cluster don't
have graphics hardware and we want to use it in parallel. See
`https://www.paraview.org/Wiki/Setting_up_a_ParaView_Server <https://www.paraview.org/Wiki/Setting_up_a_ParaView_Server>`__
for more information.

The script can be submitted via

.. code-block:: bash

   llsubmit start_paraview.sh

You then have to wait a moment and check for the pvservers to be ready
to connect.

.. code-block:: bash

   tail -f job.srv23ib.831086.out

Please replace the job name by the one you got when you submitted the
job.

If you can see something like

.. code-block:: bash

   Waiting for client...
   Connection URL: cs://i20r01c02s09:11111
   Accepting connection(s): i20r01c02s09:11111

your pvserver is up and running and ready to connect. Check your running
job with ``llu`` to get the complete ID of the leading compute node your
job is running on, e.g. ``i20r01c02s09ib``.

Setup remote visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you have already used remote visualization of LRZ you can find a file
``.vis_job.ll`` in your home directory. Open it and modify it to use
Intel MPI, i.e. set ``#@ job_type = MPICH`` instead of
``#@ job_type = parallel``. It should look like this:

.. code-block:: bash

   #!/bin/bash
   #####  job_type = parallel
   #@ job_type = MPICH
   #@ class = vis
   #@ node = 1
   #@ tasks_per_node=28
   #@ wall_clock_limit = 2:00:00
   #@ network.MPI = sn_all,not_shared,us
   #@ notification=never
   #@ node_usage = not_shared
   #@ island_count=1,1
   #@ node_topology=island
   #@ initialdir = .
   #@ output = vncjob.$(schedd_host).$(jobid).out
   #@ error = vncjob.$(schedd_host).$(jobid).err
   #@ energy_policy_tag = testtag
   #@ minimize_time_to_solution = yes
   #@ notification=never
   #@ node_resources = ConsumableMemory(16GB)
   #@ queue
   . /etc/profile . /etc/profile.d/modules.sh
   hostname
   /opt/TurboVNC/bin/vncserver -geometry 1280x900
   sleep 48h

Submit the job with ``llsubmit .vis_job.ll``.

Use ``cat vncjob.srv23ib.831117.err`` and look for something like this:

.. code-block:: bash

   Desktop 'TurboVNC: vis01:2 (username)' started on display vis01:2

This tells you which node and which display you have to use for
connecting with your VNC viewer. Start the VNC viewer with

.. code-block:: bash

   vncviewer -via username@hw.supermuc.lrz.de vis01:2

Now you have a nice GUI on the visualization node. Open a Terminal and
load the right modules:

.. code-block:: bash

   module rm poe mpi.ibm
   module load mpi.intel paraview/5.2.0
   unset I_MPI_DEVICE

Start the ParaView client with ``vglrun paraview``. Klick on ``connect``
and enter a new server. The host must be the leading compute node from
above, in this example, it is ``i20r01c02s09ib``. The port is ``11111``.
When you hit the connect button in the menu, you should have access to
all the resources you asked for in your job script and are ready to open
your data.

keyboard issue using vncviewer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A common problem is that the keyboard mapping gets all mixed-up after
vncviewer windows is deselected. To avoid this problem, add in
~/.vnc/xstart before running vncviewer:

.. code-block:: bash

   export XKL_XMODMAP_DISABLE=1

