Parameter File
==============

A parameter file for SeisSol with all possible input parameters, its
default values and a short description of these values can be found in
the [parameter file ]
(`https://github.com/SeisSol/parameter-file/blob/master/parameters.par <https://github.com/SeisSol/parameter-file/blob/master/parameters.par>`__).
For detailed information about each section, please see below.

General Information
-------------------

The parameter file in SeisSol is based on the Fortran NAMELIST format.
The file is divided into the different sections described below. Each
section has a set of configuration parameters that influences the
execution of SeisSol. Each configuration parameter can be one of the
following:

-  **Integer**
-  **Float**
-  **Array** Arrays can contain integers or floats and have a fixed
   size. The elements are separated by spaces (``1 2 3 4``)
-  **String**
-  **Path** A path references a file or a directory and is given as a
   string in the parameter file with additional restrictions. Paths
   might be an absolute or relative to the starting directory of the
   execution. If the path is used for an output file, the user has to
   make sure that the directory exists. (E.g. If the path is set to
   "output/wavefield", then the directory "output" must exist.)
-  **Path prefix** Path prefixes are similar to paths. However, SeisSol
   will automatically append a filename extension or other suffixes
   (e.g. "-fault").

Sections
--------

Additional, more detailed information on several section are listed
here.

Equations
~~~~~~~~~

IniCondition
~~~~~~~~~~~~

Boundaries
~~~~~~~~~~

DynamicRupture
~~~~~~~~~~~~~~

FL
^^

Friction constitutive law. 2: linear slip-weakening

Reference point
^^^^^^^^^^^^^^^

The slip rate is defined as the velocity difference of the two sides of
a fault, that is,

:math:`\Delta v=v^{+}-v^{-}`.

A practical issue is to define which side at an interface corresponds to
"+" and which one to "-". The reference point defines which side is
which and it is **crucial** to set it correctly.

The parameters ``XRef, YRef, ZRef`` define the coordinate vector of the
reference point, which we denote with **r**. Furthermore, the
``refPointMethod`` has to specified, whose effect is outlined in the
following.

1. | **``refPointMethod=0``**
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
   | This method works, as long as sign of the first dot product is the
     same for all faces tagged as being part of the fault.
   | *Example:* One has a planar fault with normal **N** and an
     arbitrary point **z** on the plane. Then a good reference point
     would be **z**\ +\ **N**. In this case, **n**\ =\ **N** and
   | :math:`\left<\mathbf{z}+\mathbf{N}-\mathbf{x},\mathbf{n}\right>=\left<\mathbf{z}-\mathbf{x},\mathbf{N}\right>+\left<\mathbf{N},\mathbf{N}\right>=\left<\mathbf{N},\mathbf{N}\right>`
   | that is, the first dot product becomes independent of the face.

2. | **``refPointMethod=1``**
   | Here, **r** should be rather be called reference direction instead
     of reference point.
   | We define, with **n** again being a face's normal,
   | :math:`\text{isPlus}:=\left<\mathbf{r},\mathbf{n}\right>>0`
   | *Example:* One has a planar fault with normal **N**. Then a good
     reference direction would be **N**.

*Application Example:* Assume you have chosen a *enu* coordinate system
(x=east, y=north, z=up). Your fault is in the x-z-plane with y=0
(strike-slip fault) and you set the reference point to (0,10000,0) with
``refPointMethod=0``. Then, the faces with normal (0,-1,0) make up the
"+"-side. In this case, all vertices of the "+"-tetrahedron lie in the
half-space :math:`y\ge 0`.

In the fault output, a strike and dip direction is defined (see
``create_fault_rotationmatrix.f90``). For the normal (0,-1,0), one would
obtain (-1,0,0) as strike direction (west). Recalling the definition of
the slip rate, a positive slip rate indicates left-lateral motion.

Read the article `Left-lateral, right-lateral, normal,
reverse <Left-lateral,-right-lateral,-normal,-reverse>`__ for more
information.

Elementwise
~~~~~~~~~~~

Pickpoint
~~~~~~~~~

SourceType
~~~~~~~~~~

SeisSol support point source (Type=50) and standard rupture format
(Type=42). FileName: input file.

MeshNML
~~~~~~~

SeisSol supports three mesh format: Gambit3D-fast, Netcdf and PUML

Discretization
~~~~~~~~~~~~~~

-  **ClusteredLTS**: An integer defining the algorithm used for
   clustered local time-stepping. ``1`` means global time-stepping,
   ``2``, ``3``, ... defines the multi-rate for the time steps of the
   clusters. (Default: 0)

SpongeLayer
~~~~~~~~~~~

Output
~~~~~~

-  **OutputRegionBounds**: A floating point array of size 6 that
   describes the region of the wave field that should be written. The
   region is specified in the form ``xmin xmax ymin ymax zmin zmax``.
   The default is ``0 0 0 0 0 0`` and means that the entire domain
   should be written.

Checkpoints
^^^^^^^^^^^

Checkpoints are also configured in the output section.

-  **checkPointInterval**: The checkpoint interval is a non-negative
   floating point number and gives the interval in *simulated time* for
   checkpoints. If the interval is set to zero, no checkpoints will be
   generated (Default: 0)
-  **checkPointFile**: This parameter is a *path prefix* to the location
   of the checkpoint(s).
-  **checkPointBackend**: The checkpoint back-end is specified via a
   string. Currently, the following backends are supported: ``posix``,
   ``hdf5``, ``mpio``, ``mpio_async``, ``sionlib``, ``none``. If
   ``none`` is specified, checkpoints are disabled. To use the HDF5,
   MPI-IO or SIONlib back-ends you need to compile SeisSol with HDF5,
   MPI or SIONlib respectively. (Default: ``none``)
   **Warning**: When using an asynchronous back-end (``mpio_async``),
   you might lose **2 \* checkPointInterval** of your computation.

You cannot explicitly specify to load a checkpoint. If the active
checkpoint back-end finds a valid checkpoint during the initialization,
it will load it automatically.

The parallel checkpoint back-ends (HDF5, MPI-IO, SIONlib) support
several tuning [[environment variables]].

**Hint:** Currently only the output of the wave field is designed to
work with checkpoints. Other outputs such as receivers and fault output
might require an additional post-processing when SeisSol is restarted
from a checkpoint.

AbortCriteria
~~~~~~~~~~~~~

EndTime = 10.0

Analysis
~~~~~~~~

Debugging
~~~~~~~~~

