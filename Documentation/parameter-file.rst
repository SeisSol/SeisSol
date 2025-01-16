..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

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
   string in the parameter file with additional restrictions. A path might be absolute or relative to the starting directory of the execution. If the path is used for an output file, the user has to make sure that the directory exists. (E.g. If the path is set to
   "output/wavefield", then the directory "output" must exist.)
-  **Path prefix** Path prefixes are similar to paths. However, SeisSol
   will automatically append a filename extension or other suffixes
   (e.g. "-fault").


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

