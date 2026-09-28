..
  SPDX-FileCopyrightText: 2018 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _off_fault_receivers:

Off fault receivers
===================

Introduction
------------

Receivers can be configured using the namelist output. Here is a
commented example:

.. code-block:: Fortran

  &Output
  receiverFormat = 'csv' ! Can be 'hdf5' or 'csv'
  pickdt = 0.01 ! Pickpoint Sampling
  pickDtType = 1 ! Pickpoint Type
  RFileName = 'receivers.dat' ! Record Points in extra file
  /

If ``receiverFormat = 'csv'`` (the default), each receiver trace is written to a separate ``.dat``
file. If ``receiverFormat = 'hdf5'``, all of them go into a single file, ``-receivers.h5``; see
`The HDF5 receiver file`_ below for what it holds.

If pickDtType = 2, the output is generated every N time steps, where N is
set by pickdt. If pickDtType = 1, output is generated every pickdt
second.

receivers.dat is an ASCII file describing the coordinates of the receivers in
the form:

::

  x1 y1 z1
  x2 y2 z2
  (...)
  xn yn zn


The receivers files contain the time-histories of the stress tensor (6 variables) and the particle velocities (3).
Currently, there is no way to write only a subset of these variables.

The variable :code:`ReceiverOutputInterval` (in the section :code:`Output` of the :ref:`parameter-file`) controls the frequency of flushing receiver time-histories. If not specified, they are written at the end of the simulation.

The HDF5 receiver file
----------------------

Everything sits under an HDF5 group named ``receivers``.

What a receiver records follows from the material of the element it sits in, so
not all of them record the same quantities. The receivers are therefore gathered
into groups that share a quantity set, and each group becomes a dataset
``group0``, ``group1``, ... of its own. A run in which every receiver records the
same quantities -- the usual case -- has exactly one of them.

A dataset is indexed by sample and by receiver:

::

  /receivers/group0        (samples, receivers)   compound

One element is a whole sample of one receiver, as a compound whose members are
the quantities: ``Time`` first, then the material quantities, then the derived
ones if :code:`ReceiverComputeRotation` or :code:`ReceiverComputeStrainRate` are
on. That is the same memory as a ``(sample, receiver, quantity)`` array of
numbers, with the quantity axis named rather than numbered, so the names and the
types come out of the file itself. All the samples of one receiver lie together,
which is the access a post-processing step usually wants.

Two attributes describe a dataset:

``Quantities``
  the quantity set it holds, in the form the grouping reads back

``NumberOfPoints``
  how many receivers it holds over all ranks

Beside the datasets, and written once, are the columns describing the receivers,
in the order the ranks contributed them:

::

  /receivers/Index         (receivers, 2)    group and row within that group
  /receivers/PointId       (receivers,)      the receiver's line in the receiver file
  /receivers/Coordinates   (receivers, 3)    where it sits

The receivers are renumbered so that every rank owns one run of each group;
``Index`` is what leads from a line of the receiver file back to the row of the
dataset that holds it.

With numpy and h5py, reading the trace of the receiver on line ``n`` of the
receiver file is therefore:

.. code-block:: python

  import h5py

  with h5py.File("output-receivers.h5") as f:
      receivers = f["receivers"]
      index = receivers["Index"][:]
      row = (receivers["PointId"][:] == n).nonzero()[0][0]
      group, column = index[row]
      trace = receivers[f"group{group}"][:, column]
      time = trace["Time"]
      v1 = trace["v1"]

Storage chunking
~~~~~~~~~~~~~~~~

How far a storage chunk reaches along the sample axis can be set with
:code:`receiversamplechunk` in the :code:`Output` section; zero, the default,
lets the writer take one write as one chunk. The chunking is settled when the
file is created while the number of samples a write carries varies as soon as
:code:`pickdt` does not divide :code:`ReceiverOutputInterval`, so a run that is
read back sample-wise rather than receiver-wise may want it given.


Rotational Output
-----------------
You can additionally choose to write the rotation of the velocity field by setting :code:`ReceiverComputeRotation=1` in the parameter file.
The rotation of the vector field is defined as :math:`\text{rot} v = \begin{pmatrix} \partial_2 v_3 - \partial_3 v_2 \\ \partial_3 v_1 - \partial_1 v_3 \\ \partial_1 v_2 - \partial_2 v_1 \\ \end{pmatrix}`.

Strain Rate Output
------------------
Furthermore, you can also output the strain rate by setting :code:`ReceiverComputeStrainRate=1`.

Placing free-surface receivers
------------------------------

Placing receivers on the free-surface requires special care when a
realistic topography is used. The procedure to move receivers exactly to
the surface is described
`here <https://github.com/SeisSol/Meshing/tree/master/place_receivers>`__.

Compiling place_receivers on SuperMUC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Load the relevant :ref:`modules <compile_run_supermuc>`.

.. code-block:: bash

  git submodule update --init
  mkdir build && cd build
  cmake ..
  make -j
