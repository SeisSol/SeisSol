..
  SPDX-FileCopyrightText: 2018 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Known issues
============

Asynchronous output
-------------------

Some MPI versions have issues (i.e. deadlocks or crashes)
with threading and IO functionality; especially the Hdf5 library.

The MPI-IO
implementation is either not threadsafe or too strict when it comes to
locking. This can lead to crashes or deadlocks when using the
asynchronous output.

A workaround is to either use the POSIX backend for the Xdmf writer module,
or to set ``ASYNC_MODE=sync`` instead of ``ASYNC_MODE=thread``.

.. _"holes"-in-the-fault-output:

"Holes" in the fault output
---------------------------

"Holes" in the fault output may appear if the reference point or
reference vector, defined by (Xref, Yref, and Zref) is wrongly located
(point exactly on the fault, or null vector).
