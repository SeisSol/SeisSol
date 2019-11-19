Point source (older implementation)
===================================

Using ``Type = 50`` is the old and non-optimal way to include a point
source in SeisSol. It might nevertheless still be useful for modeling a
non double-couple point source (not currently possible with the nrf
source description).

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

Using this implementation one can specify point sources in the following form:

.. math ::
  \begin{aligned}
  \frac{\partial}{\partial t}\sigma_{ij} + \dots &= M_{ij} \cdot S_k(t)\cdot \delta(x - \xi_k) \\
  \rho \frac{\partial}{\partial t} u_i + \dots &= d_i \cdot S_k(t) \cdot \delta(x - \xi_k)
  \end{aligned}

For details about the source term read section 3.3 of `An arbitrary high-order discontinuous Galerkin method for elastic
waves on unstructured meshes – I. The two-dimensional isotropic case with external source terms 
<https://academic.oup.com/gji/article-lookup/doi/10.1111/j.1365-246X.2006.03051.x>`__

Mij is defined in a fault local coordinate system defined by strike, dip and rake, see for instance here:
`https://github.com/SeisSol/SeisSol/blob/master/src/SourceTerm/PointSource.cpp#L48 <https://github.com/SeisSol/SeisSol/blob/master/src/SourceTerm/PointSource.cpp#L48>`__
