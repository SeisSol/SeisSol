..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Known issues
============

Download
--------

-  **The submodule directories are empty.** To clone the code including
   all submodules, use

   ::

      git clone --recursive https://github.com/SeisSol/SeisSol.git

   If you already cloned the repository without ``--recursive``, you can
   use the following command to get the submodules.

   ::

      git submodule update --init --recursive

Installation
------------

**The build process fails because of missing files.** SeisSol
requires some external libraries which are integrated via git submodules. For some of the libraries it is important to have the
correct version (e.g. the XMDF Writer) because of changes in the
interface of the library. If you clone the submodules with git, you
should always get the correct version. However, downloading the
SeisSol and its submodules directly from the Github homepage might result in an incorrect combination of versions.

Asynchronous output
-------------------

Some MPI versions have issues with threads and I/O. The MPI-IO
implementation is either not threadsafe or too strict when it comes to
locking. This can lead to crashes or deadlocks when using the
asynchronous output. A simple workaround is to use the POSIX back-end of
the XDMF writer, this is compiling SeisSol without HDF5. For
checkpointing, one should also switch to POSIX or SIONlib.

.. _easi-and-intel/16.0:

easi and Intel/16.0
-------------------

Using Intel/16.0 in a local cluster, compile error occurs in
PUMLreader.cpp Using Intel/17.0 will solve this.

Using Intel/16.0 on SuperMUC, we observed some bugs in the fault stress
initialization (you can see the imprint of some partitions in the
initial stress). The bugs are not showing off with latest Intel/17.0
module.

.. _"holes"-in-the-fault-output:

"Holes" in the fault output
---------------------------

"Holes" in the fault output may appear if the reference point or
reference vector, defined by (Xref, Yref, and Zref) is wrongly located
(point exactly on the fault, or null vector).
