..
  SPDX-FileCopyrightText: 2018 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

.. _parameter-file:

Parameter file
==============

General information
-------------------

The parameter file in SeisSol is based on the Fortran NAMELIST format.
The file is divided into different sections. Each section has a set of configuration parameters that influences the
execution of SeisSol. Each configuration parameter can be one of the
following:

-  **Integer**
-  **Float**
-  **Array** Arrays can contain integers or floats and have a fixed
   size. The elements are separated by spaces (``1 2 3 4``)
-  **String**
-  **Path** A path references a file or a directory and is given as a
   string in the parameter file with additional restrictions. A path might be absolute or relative to the starting directory of the execution.
-  **Path prefix** Path prefixes are similar to paths. However, SeisSol
   will automatically append a filename extension or other suffixes
   (e.g. "-fault"). SeisSol creates the directory of the output prefix
   (``OutputFile``), with its parents, if it does not exist yet; e.g. for
   "runs/a/output/data", it creates "runs/a/output".


.. _parameter-overrides:

Overriding parameters on the command line
-----------------------------------------

``--set section.field=value`` gives a parameter on the command line in place of the
value in the parameter file; it can be given several times, and a later one wins
over an earlier one. Section and field are written as in the parameter file, in any
case, and the value is read as if it stood there, quotes and Fortran exponents
included:

.. code-block:: bash

    ./SeisSol_YOUR_CONFIGURATION parameters.par \
        --set output.outputfile=runs/strong-16/data \
        --set abortcriteria.endtime=2.5d0 \
        --set output.wavefieldoutput=0

That way, a series of runs can share one parameter file. A section the file does not
have is added. An override of a field SeisSol does not read is an error, unlike a
field in the file that SeisSol does not know, which only gives a warning: a
misspelt override would otherwise silently leave the run unchanged.

.. _parameter-check:

Checking a run before starting it
---------------------------------

``--check`` reads the parameter file and the overrides as a run does, but instead of
running, it checks what a run finds out only later, when it may have waited in a queue
for hours:

* that the files the parameters name exist, those that the run opens (the receiver
  file, e.g., only if the receiver output is on);
* that the YAML files among them, i.e. the easi models, parse, though not what they
  hold, nor the files they refer to in turn;
* that the mesh opens, and how many cells and vertices it has;
* that the output directory can be written, or created;
* that the checkpoint given with ``-c`` exists.

Then it reads the whole mesh once, cell by cell, and checks what the run would check
only after reading and partitioning it:

* that the receivers lie in cells of the mesh, the point sources as well, and that the
  on-fault receivers lie on the fault, i.e. across a fault face along its normal; it
  takes the points with the same tests as the run. For a point that lies in no cell, it
  tells whether the point lies above the free surface, and by how much, or outside of
  the mesh. It warns of an on-fault receiver that the run moves onto the fault from
  farther away than the fault face is large;
* that the face map (``facemapfile``, or the default one) knows the boundary condition
  tags of all faces, and that the build supports the boundary conditions that the cells
  of each configuration have;
* that the cells are oriented as the run requires, which a cell without volume, or a
  scaling of the mesh that mirrors it, is not.

It warns, too, if dynamic rupture is on but the mesh has no fault faces, if the
free-surface output is on but the mesh has no faces of the free surface, and if the
configuration map names a group that no cell of the mesh has.

It runs on one process, without MPI and without devices, with the threads of OpenMP,
e.g. on a login node with the executable built for the compute nodes:

.. code-block:: bash

    ./SeisSol_YOUR_CONFIGURATION --check=check.json parameters.par \
        --set output.outputfile=runs/strong-16/data

Reading the whole mesh holds its vertices in memory, 24 bytes each, and a chunk of its
cells; for a mesh of 10 million cells, it took 2 s on two cores. ``--check-depth=files``
leaves it out, ``--check-depth=mesh`` is the default; ``--check-depth`` alone implies
``--check``.

The login node has to be able to start the executable. One built with CUDA starts
without the NVIDIA driver as well; where an executable does not start, ``ldd`` names the
libraries the node lacks, e.g. those that a GPU-aware MPI library needs. An executable
built for the instructions of the CPUs of the compute nodes may stop on a login node
with other CPUs.

An error in the parameters stops the check as it would stop the run; what the check
itself finds, it collects: as problems, which keep the run from starting or from going
as intended, and as warnings. The exit code is 0 if the run can start, 1 otherwise.
``--check=<file>`` writes a report as JSON, also when an error or a crash stops the
check: its ``status`` is ``passed`` or ``failed``, or ``terminated`` if a request to
stop reached it (then ``terminatedBy`` gives the signal); ``problems`` and ``warnings``
list what the check found; ``stoppedAt`` names the step an error or a crash stopped
(``parameters``, ``files``, ``yaml``, ``mesh``, ``checkpoint``, ``output``, ``points``
or ``meshScan``), and ``error`` gives its ``message`` and its ``time``, for an error its
``source`` in the code (``file``, ``line``, ``function``), for a crash the ``signal``
(see :ref:`run-summary`); ``depth`` is the depth of the check; ``mesh`` has the
``cells`` and ``vertices`` of a PUML mesh, and with the whole mesh read, ``scan``: its
bounding box (``min``, ``max``), the faces by their boundary condition tag (``tags``)
and face type (``faceTypes``), the cells by their group (``groups``), the number of
fault faces (``faultFaces``), the radii of the inspheres of the smallest and the largest
cell (``smallestInsphere`` of the cell ``smallestCell`` of the file, and
``largestInsphere``), and how long reading it took (``seconds``) with how much memory
(``bytes``); ``points`` has the ``receivers``, ``faultReceivers`` and ``sources``, each
with the ``file``, the ``count`` of points, how many were ``found``, and the ``missing``
ones (their ``index``, counted from zero as the run counts them, the ``point``, the
height of the free surface above or below it, ``surface``, and a ``hint``), and for the
on-fault receivers the ones ``moved`` farther than their face is large; ``files`` lists
the files the parameters name, whether they exist and whether the run opens them;
``output`` gives the output directory, whether it is ``writable``, and which outputs are
``enabled``; the report has the ``build``, the ``invocation``, the ``configurations``,
the ``endTime`` and the ``unknownFields`` of the parameter file as well.

Commented parameter file
------------------------

.. literalinclude:: parameters.par
   :language: fortran

Sections
--------

Additional, more detailed information on several sections are listed
here.

Dynamic rupture
~~~~~~~~~~~~~~~

Reference point
^^^^^^^^^^^^^^^

The slip rate is defined as the velocity difference between the two sides of
a fault, that is,

:math:`\Delta v=v^{+}-v^{-}`.

A practical issue is to define which side at an interface corresponds to
"+" and which one to "-". The reference point defines which side is
which and it is **crucial** to set it correctly.

The parameters ``XRef, YRef, ZRef`` define the coordinate vector of the
reference point, which we denote with **r**. Furthermore, the
``refPointMethod`` has to specified, whose effect is outlined in the
following.

1. | **refPointMethod=0**
   | In order to decide if a side of a fault is "+" or "-" we compute
     the vectors **n**, **x**, and **y** for a face of a tetrahedron,
     whose boundary condition indicates it to be part of the fault. The
     vector **n** is the face's normal, which always points outward with
     respect to the tetrahedron. The vector **x** is an arbitrary vertex
     of the face and the vector **y** is face's the missing vertex, that
     is, the vertex which belongs to the tetrahedron but not to the
     face.
   | We define
   | :math:`\text{isPlus}:=\left<\mathbf{r}-\mathbf{x},\mathbf{n}\right>\cdot\left<\mathbf{y}-\mathbf{x},\mathbf{n}\right>>0`
   | isPlus is only true whenever **r**-**x** and **y**-**x** point in
     the same direction (lie in the same half-space w.r.t. **n**).
   | This method works, as long as the sign of the first dot product is the
     same for all faces tagged as being part of the fault.
   | *Example:* One has a planar fault with normal **N** and an
     arbitrary point **z** on the plane. Then a good reference point
     would be **z**\ +\ **N**. In this case, **n**\ =\ **N** and
   | :math:`\left<\mathbf{z}+\mathbf{N}-\mathbf{x},\mathbf{n}\right>=\left<\mathbf{z}-\mathbf{x},\mathbf{N}\right>+\left<\mathbf{N},\mathbf{N}\right>=\left<\mathbf{N},\mathbf{N}\right>`
   | that is, the first dot product becomes independent of the face.

2. | **refPointMethod=1**
   | Here, **r** should be rather be called reference direction instead
     of reference point.
   | We define, with **n** again being a face's normal,
   | :math:`\text{isPlus}:=\left<\mathbf{r},\mathbf{n}\right>>0`
   | *Example:* One has a planar fault with normal **N**. Then a good
     reference direction would be **N**.

*Application Example:* Assume you have chosen a *enu* coordinate system
(x=east, y=north, z=up). Your fault is in the x-z-plane with y=0
(strike-slip fault) and you set the reference point to (0,10000,0) with
``refPointMethod=0``. Then, the faces with normal (0,+1,0) make up the
"+"-side. In this case, all vertices of the "+"-tetrahedron lie in the
half-space :math:`y\ge 0`.

In the fault output, a strike and dip direction is defined (see
``create_fault_rotationmatrix.f90``). For the normal (0,-1,0), one would
obtain (-1,0,0) as strike direction (west). Recalling the definition of
the slip rate, a positive slip rate indicates left-lateral motion.

Read the article `Left-lateral, right-lateral, normal,
reverse <Left-lateral,-right-lateral,-normal,-reverse>`__ for more
information.
