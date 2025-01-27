..
  SPDX-FileCopyrightText: 2018-2024 SeisSol Group

  SPDX-License-Identifier: BSD-3-Clause
  SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/

  SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

Point source (older implementation)
===================================

Using ``Type = 50`` is an older and less-optimized way to include a point
source in SeisSol. It might nevertheless still be useful for modeling a
non-double-couple point source, which is not currently possible with the nrf
source description.

Add the following section to your parameter file:

::

   &SourceType
   Type = 50
   FileName = 'source.dat'
   /

Where source.dat has the following format

::

   header line (e.g. 'Seismic Moment Tensor')
   M11 M12 M13
   M21 M22 M23
   M31 M32 M33
   header line (optional, but has to contain 'velocity' to be recognized)
   d1 d2 d3    (optional)
   header line (e.g. 'Number of subfaults')
   nsubfaults
   header line (e.g. 'x y z strike dip rake area Onset time')
   x y z strike dip rake area Onset time(1)
   x y z strike dip rake area Onset time(2)
   ....
   x y z strike dip rake area Onset time(nsubfaults)
   header line (e.g. 'source time function')
   dt ndt
   header line (e.g. 'samples')
   STF(1,1)
   STF(1,1)
   ...
   STF(1,ndt)
   STF(2,1)
   ...
   STF(nsubfaults,1)
   ...
   STF(nsubfault,ndt)

This format describes point sources following the equations

.. math ::
  \begin{aligned}
  \frac{\partial}{\partial t}\sigma_{xx} - (\lambda + 2\mu) \frac{\partial}{\partial x} u - \lambda \frac{\partial}{\partial y} v - \lambda \frac{\partial}{\partial z} w &= M_{xx} \cdot S_k(t)\cdot \delta(x - \xi_k) \\
  \frac{\partial}{\partial t}\sigma_{yy} - \lambda \frac{\partial}{\partial x} u - (\lambda+2\mu) \frac{\partial}{\partial y} v - \lambda \frac{\partial}{\partial z} w &= M_{yy} \cdot S_k(t)\cdot \delta(x - \xi_k) \\
  \frac{\partial}{\partial t}\sigma_{zz} - \lambda \frac{\partial}{\partial x} u - \lambda \frac{\partial}{\partial y} v - (\lambda + 2\mu)  \frac{\partial}{\partial z} w &= M_{zz} \cdot S_k(t)\cdot \delta(x - \xi_k) \\
  \frac{\partial}{\partial t}\sigma_{xy} - \mu \frac{\partial}{\partial x} v - \mu \frac{\partial}{\partial y} u &= M_{xy} \cdot S_k(t)\cdot \delta(x - \xi_k) \\
  \frac{\partial}{\partial t}\sigma_{yz} - \mu \frac{\partial}{\partial z} v - \mu \frac{\partial}{\partial y} w &= M_{yz} \cdot S_k(t)\cdot \delta(x - \xi_k) \\
  \frac{\partial}{\partial t}\sigma_{xz} - \mu \frac{\partial}{\partial z} u - \mu \frac{\partial}{\partial x} w &= M_{xz} \cdot S_k(t)\cdot \delta(x - \xi_k) \\
  \rho \frac{\partial}{\partial t} u - \frac{\partial}{\partial x} \sigma_{xx} - \frac{\partial}{\partial y} \sigma_{xy} - \frac{\partial}{\partial z} \sigma_{xz} &= d_x \cdot S_k(t)\cdot \delta(x - \xi_k) \\
  \rho \frac{\partial}{\partial t} v - \frac{\partial}{\partial x} \sigma_{xy} - \frac{\partial}{\partial y} \sigma_{yy} - \frac{\partial}{\partial z} \sigma_{yz} &= d_y \cdot S_k(t)\cdot \delta(x - \xi_k) \\
  \rho \frac{\partial}{\partial t} w - \frac{\partial}{\partial x} \sigma_{xz} - \frac{\partial}{\partial y} \sigma_{yz} - \frac{\partial}{\partial z} \sigma_{zz} &= d_z \cdot S_k(t)\cdot \delta(x - \xi_k). \\
  \end{aligned}

For details about the source term, we refer to Section 3.3 of `An arbitrary high-order discontinuous Galerkin method for elastic
waves on unstructured meshes – I. The two-dimensional isotropic case with external source terms
<https://academic.oup.com/gji/article-lookup/doi/10.1111/j.1365-246X.2006.03051.x>`__

Mij and di are defined in a fault local coordinate system defined by strike, dip and rake, see for instance here:
`https://github.com/SeisSol/SeisSol/blob/master/src/SourceTerm/PointSource.cpp#L48 <https://github.com/SeisSol/SeisSol/blob/master/src/SourceTerm/PointSource.cpp#L48>`__

The above equations also hold for viscoelastic or anisotropic materials.
In the viscoelastic case, the equations are extended by the memory variables.
In the anisotropic case, :math:`\lambda` and :math:`\mu` are replaced by the entries of the Hooke tensor :math:`c_{ij}`.

For poroelastic materials, we add the possibility to consider forces in the fluid or pressure sources.
To do so, add these two lines before ``Number of subfaults``:

::

   header line (optional, but has to contain 'pressure' to be recognized)
   p
   header line (optional, but has to contain 'fluid' to be recognized)
   f1 f2 f3

In the case of a poroelastic material, the equations of motion differ from above equations. For the velocities in :math:`x` direction, they read:

.. math ::
  \begin{aligned}
   \frac{\partial}{\partial t} u - \frac{1}{\rho^{(1)}} \frac{\partial}{\partial x} \sigma_{xx} - \frac{1}{\rho^{(1)}} \frac{\partial}{\partial y} \sigma_{xy} -  \frac{1}{\rho^{(1)}} \frac{\partial}{\partial z} \sigma_{xz} - \frac{1}{\rho^{(1)}} \frac{\rho_f}{m}\frac{\nu}{\kappa} u_f &= d_x \cdot S_k(t) \delta(x - \xi_k) \\
   \frac{\partial}{\partial t} u_f - \frac{1}{\rho^{(2)}} \frac{\partial}{\partial x} \sigma_{xx} - \frac{1}{\rho^{(2)}} \frac{\partial}{\partial y} \sigma_{xy} -  \frac{1}{\rho^{(2)}} \frac{\partial}{\partial z} \sigma_{xz} - \frac{1}{\rho^{(2)}} \frac{\rho}{\rho_f}\frac{\nu}{\kappa} u_f &= f_x \cdot S_k(t) \delta(x - \xi_k) \\
  \end{aligned}

with

.. math ::
  \rho^{(1)} &= \left(\rho - \frac{\rho_f^2}{m}\right), \quad \rho^{(2)} &= \left(\rho_f - \frac{m \rho}{\rho_f} \right).\\

In particular, the terms :math:`d_x` and :math:`f_x` are not divided by :math:`\rho^{(1)}` or :math:`\rho^{(2)}`, respectively.
