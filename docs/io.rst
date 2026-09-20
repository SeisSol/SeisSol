..
  SPDX-FileCopyrightText: 2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _io_infrastructure:

The IO Infrastructure
=====================

SeisSol uses a single IO module for everything it writes. It describes a write
rather than performing it, hands the description to a writer that may run on a
thread or a process of its own, and carries on computing. Everything a
simulation produces goes through it:

* the mesh outputs: the wavefield, the free surface and the elementwise fault
  output, at degree zero or as high-order VTKHDF
* the point outputs: the off-fault receivers and the on-fault receivers
* checkpointing of arbitrary friction laws and equations
* the tables a run leaves beside its output, such as the clustering and the
  thread pinning

How a write is described
------------------------

A write is a **plan**: a list of instructions, each naming a place in a file, a
name to write under, and a **data source** to take the values from. The plan is
built while the simulation holds still at a synchronization point, and it is the
plan -- not the data -- that travels to the writer.

A data source is either carried inside the plan, for the few values that are
small enough, or it is a buffer the writer copies out before the simulation
resumes. Which of the two it is follows from where the memory comes from, and
nothing else in the module depends on it.

Shapes
~~~~~~

Every data source carries its **shape**: what its dimensions are, and what each
of them means to the file. Two properties are asked of each dimension, and they
are independent of each other.

How the ranks share it:

``Replicated``
  Every rank holds the same entries, and they are written once.

``Distributed``
  The entries are split across the ranks and concatenated in rank order. Its
  size is not declared -- how many entries a rank contributes is only known once
  the data is there, and the total follows from a scan over the ranks.

How it grows from one write to the next:

``Fixed``
  The dataset has this extent once and for all. Writing the same dataset twice
  is an error.

``Appended``
  The dataset is unlimited here, and every write adds its entries at the end. A
  dimension carries the extent of *one* write, so how far the dataset grows, not
  how large it ends up.

At most one dimension of a shape is distributed and at most one is appended, but
neither has to be the first one, and they may be the same one. Two shapes worth
knowing, because they are what the outputs use:

* A VTKHDF time series concatenates its steps along the same dimension the ranks
  are split along -- one flat array that the step offsets slice. That dimension
  is distributed *and* appended.
* A table of samples per point grows along one dimension and is distributed
  along another, with the quantities of a sample fixed in a third.

A dataset grows by what its appended dimension carries: the rank total where
that dimension is also distributed, and its declared size otherwise. How far a
storage chunk reaches can be given per dimension; by default one write is one
chunk along an appended dimension, and a few megabytes along a distributed one.
The chunking is settled when the dataset is created and holds for the whole run,
so an output whose writes carry a varying number of entries wants it given
rather than derived from whichever number the first write happened to have.

Mesh output
-----------

The wavefield, the free surface and the elementwise fault output share one
writer. It is fed the points of a cell and the values on them, and writes either
VTKHDF -- when a ``vtkorder`` is set for that output -- or Xdmf.

At degree zero the corners a cell shares with its neighbours are written once
and the cells index into them, which makes the point array as large as the mesh
rather than as large as the mesh times the number of cells a vertex touches. The
merging happens within a rank, so points on a partition boundary stay
duplicated, which is what a reader going partition by partition expects. It can
be switched off, see :ref:`SEISSOL_IO_VERTEXFILTER <env_vars>`. From degree one
on the points are Lagrange nodes, which neighbouring cells deliberately do not
share, and nothing is merged.

File groupings
~~~~~~~~~~~~~~

Every mesh output writes one self-contained file per output step by default. Two
other groupings are available, and both need the output to be written as VTKHDF:

``snapshot``
  One file per step, each complete in itself. The default.

``incremental``
  The data that does not change -- the points and the connectivity -- is written
  once per run, into a file of its own, and every snapshot references it. Worth
  it whenever the mesh is large against the fields written on it.

``monolith``
  One file for the whole run, as a VTKHDF time series. The step bookkeeping
  lives inside that file, so no ``.pvd`` collection is written alongside it.

``outputtimeseries`` in the ``output`` section picks one for all three mesh
outputs at once. ``wavefieldtimeseries`` and ``surfacetimeseries`` in the same
section, and ``timeseries`` in the ``elementwise`` section, override it for one
of them. An output written as Xdmf says so and falls back to ``snapshot``,
rather than ending the run over a setting that a different parameter --
``vtkorder`` -- decides.

Point output
------------

The receivers, on-fault and off-fault, are written as tables. What a point
records follows from the material of the element it sits in, so not all of them
record the same quantities; the points are therefore gathered into groups that
share a quantity set, and every group becomes a dataset of its own rather than
one wide table with holes in it. Which groups exist is agreed across the ranks,
since a rank takes part in declaring a dataset even when it holds no point of
that group.

A dataset is indexed by sample and by point, and one element holds a whole
sample of one point as a compound of its quantities. That is the same memory as
a ``(sample, point, quantity)`` array of numbers, with the quantity axis named
instead of numbered: a reader takes the names and the types out of the file, and
a group whose quantities are not all of one type still fits. The sample axis
grows with every write and the point axis is split across the ranks, so all the
samples of one point lie together -- which is how these files are read
afterwards, one point at a time.

Alongside the tables, and written once, there is a point map giving for every
point the group it went into and its row in that group, plus whatever the output
knows about a point: its number in the receiver file, its coordinates, and for
the on-fault receivers the face and the cells on either side of it. The points
are renumbered so that every rank owns one run of each group, which is what
makes the layout expressible as a single distributed dimension; without the map
there would be no way back from a row to the point it belongs to.

See :ref:`off_fault_receivers` and :ref:`fault_receivers` for how to enable
these and what the files hold.

Tables beside the output
------------------------

The smaller tables a run leaves next to its output -- the clustering, the thread
pinning, the load balance of the initial measurement -- are written as CSV
through the same table writer. The columns are declared with their types, so the
header follows from them and the escaping of values holding a delimiter or a
quote is done in one place. These are assembled by one rank and written
directly, since nothing in the simulation waits for them.

Tuning
------

Which settings pay off depends entirely on the file system, so there is nothing
sensible for SeisSol to default to. The knobs are environment variables, and
they are listed under :ref:`env_vars`: the alignment of the bulk data, the block
size, and the MPI-IO hints passed through to the implementation.

ASYNC
~~~~~

ASYNC is the backbone of the IO infrastructure, and enables the offloading of IO
operations to dedicated threads or even processes. Internally, the IO system
serializes all write requests (to binary or Hdf5 files) and sends them to the
ASYNC executors. See :ref:`asynchronous-output` for more information.
