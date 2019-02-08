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

From the source, it seems that Mij is defined in a fault local
coordinate system defined by strike, dip and rake, see for instance
here:
`https://github.com/SeisSol/SeisSol/blob/master/src/SourceTerm/PointSource.cpp#L45 <https://github.com/SeisSol/SeisSol/blob/master/src/SourceTerm/PointSource.cpp#L45>`__
