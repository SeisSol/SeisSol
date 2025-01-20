..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Left lateral, right lateral, normal, reverse
============================================

Assume we are given a reference point **r** and a normal **n** pointing
from the "+"-side to the "-"-side of a fault segment. How does one
describe the motion of the fault?

Left-lateral and right-lateral
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"In a left-lateral (right-lateral) fault, an observer on one of the
walls will see the other wall moving to the left (right)." [J. Pujol,
Elastic Wave Propagation and Generation in Seismology]

Assume we stand on the "-"-side and look towards the "+"-side, then if
the strike vector **s** points to left, we have a left-lateral motion
(for a positive slip-rate). We formalize "points to the left" with:

:math:`l:=u\times(-n)`

where **u** is the unit vector which points up (e.g. (0,0,1) for *enu*
or (0,0,-1) for *ned*).

In SeisSol, the strike vector is (not normalized)

:math:`s:=(-e_3)\times n`

So, e.g., for *enu* we always have a left-lateral motion, as **s** and
**l** are parallel, and for *ned* we always have a right-lateral motion
as **s** and **l** are anti-parallel.

Normal and reverse
~~~~~~~~~~~~~~~~~~

"The foot wall (hanging wall) is defined as the block below (above) the
fault plane. (...) the hanging wall moves up with respect to the foot
wall and the fault is known as *reverse*. (...) the opposite happens and
the fault is said to be *normal*." [J. Pujol, Elastic Wave Propagation
and Generation in Seismology]

In SeisSol, the dip vector is (not normalized)


:math:`d:=n\times s=n\times(-e_3\times n)=-e_3+n_zn`

We used Grassmann's Identity for the last step. In particular, we
observe that the dip vector **d** is independent of the reference point,
as we obtain the same vector for -**n** and, as **n** is normalized and

:math:`d_z:=-e_3+n_z^2`,

the dip vector always points in -z direction. That is, the "+"-side
moves down for *enu* and the "+"-side moves up for *ned* (assuming
positive slip rate).

Normal or reverse depends also on the reference point. If the reference
point is inside the hanging wall (foot wall) the "+"-side corresponds to
the hanging wall (foot wall).

In summary, we obtain the following table

======== ============= ================
\        foot wall= +  hanging wall = +
======== ============= ================
z = up   reverse       normal
z = down normal        reverse
======== ============= ================

or logically

:math:`\text{isNormal}:=(+=\text{hanging wall})\leftrightarrow(z=\text{up})`

Example
~~~~~~~

We have a 60° dipping normal fault with 90° strike (West-East direction, pointing
towards the East with *enu* convention) and 0° rake. The normal of the fault
plane, which points from the foot wall to the hanging wall, is given by

:math:`n:=\frac{1}{2}\begin{pmatrix}0 & -\sqrt{3} & 1\end{pmatrix}`

Such normal vector is obtained as the result of the cross product between the dip vector

:math:`d:=\frac{1}{2}\begin{pmatrix}0 & -1 & -\sqrt{3}\end{pmatrix}`

and the respective strike vector

:math:`s:=\begin{pmatrix}1 & 0 & 0\end{pmatrix}`

Hence, we set the reference point to **x** + a **N**, where a > 0 and
**x** is an arbitrary point on the fault. In this case, the reference
point is inside the hanging wall and we obtain a normal fault.

Warning
~~~~~~~

The rake angle describes the direction of slip. In SeisSol, the convention for the rake angle is to assume a positive rake angle refers to right-lateral strike slip.
!Warning!
This is opposite to the common convention in seismology which assumes positive rake angle implying left-lateral strike slip faulting. Thus, e.g., the case of Ts0>0 and Td0=0 means initial shear stress loading in rake 180° direction.
