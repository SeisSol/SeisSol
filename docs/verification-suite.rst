..
  SPDX-FileCopyrightText: 2026 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

The verification suite
======================

``tests/verification`` runs SeisSol from end to end on small meshes that it
generates itself (``minimesh.py``), in each configuration of the executable
that ``SEISSOL_VERIFICATION_CONFIGS`` names (see :doc:`build-parameters`), and
holds the runs against exact solutions, against each other and against stored
references. This page describes how to run the suite and which runs it holds
against each other; :doc:`verification-exact-solutions` describes the cases of
dynamic rupture.

Running the suite
-----------------

The cases are tests of a build with ``TESTING=ON``, named after their
configuration, as ``verification.smoke-box48<elastic-o6-f64>``, with labels by
what they hold: ``ctest -L verification`` runs all of them, ``ctest -L dr``,
``-L exact`` or ``-L identity`` a part of them, and
``ctest -LE spectral`` leaves out the cases that take minutes each.

Identities
----------

Some runs have to agree with others, up to rounding at most, as they pose the
same problem in another form:

``identity-fault48-off``, ``identity-fault48-nofault``
   A fault tag without a fault file is no fault: its faces are faces as any
   other, bit for bit. With the friction law ``FL=0`` and a fault file
   (``fault-welded.yaml``), the fault is welded, which is the absence of a
   fault as well, up to the rounding of the flux across it.

``identity-*-relabel``
   SeisSol numbers the vertices of a cell by their ids, which decide which
   face of each of its two cells a face is, and with it the variant of every
   kernel that works on the face. Of the 14 pairs of faces that this numbering
   allows on the two sides of a face, a mesh numbered along its lattice holds
   two; its copy ``*-relabel`` numbers the vertices at random, and holds up to
   all 14 (``RELABELLED`` in ``minimesh.py``). A welded fault on it is the
   absence of a fault, and with a fully symmetric quadrature rule on the faces,
   as Dunavant's, a rupture on it is the rupture on the original, up to
   rounding.

``identity-*-mirror``
   The reference point of the fault (``XRef``) on its other side turns its
   normal around, and with it its strike, but not its dip
   (:doc:`left-lateral-right-lateral-normal-reverse`). The jump across the
   fault turns around as well, so that the fault receivers report the same
   along the strike and across the fault, and the opposite along the dip and
   the normal velocity.

``identity-rupture-*-energy-*``
   The energy output sums over the cells and the faces of a rank and then
   over the ranks: on another number of ranks or threads, in another order.

``bitwise-*-restart``
   A run continued from a checkpoint is the run without it, bit for bit; with
   rate and state and thermal pressurization, the checkpoint carries the state
   of the diffusion across the shear zone as well.

The differences of the receivers are relative to the largest value of their
column, as a value that passes through zero, as the shear stress of a wave
does, would magnify a difference by rounding beyond any bound.
