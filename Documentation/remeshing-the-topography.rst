Remeshing the topography
========================

Nowadays, high resolution topographic and bathymetric data are most of
the time available. Processing this large amount of data can really be a
challenge. For example in Gocad, intersecting such surface with other
surfaces can be time consuming and error prone. To overcome these kind
of difficulties, an idea is to coarsen the meshed topography where a
fine resolution is not necessary, before working on the surfaces in
Gocad.

As an illustration, we will process a netcdf file from GEBCO
(`http://www.gebco.net/ <http://www.gebco.net/>`__). It has been used in
Sumatra related simulations. It features a 400m resolution regular grid.
Using the script createGOCADTSurfNXNY_netcdf.py, available
`here <https://github.com/SeisSol/Meshing/tree/master/GocadRelatedScripts>`__,
we will downsample the data overall, project them and isolate a square
region away from which a fine discretisation is not necessary.

| ``python createGOCADTSurfNXNY_netcdf.py data/GEBCO_2014_2D_90.0_1.5_97.0_14.5.nc  trash.stl --subsample 2 --proj EPSG:32646  --hole 94 95 8 10``
| Now we can import the data in SimModeler5 (the version is important,
  as SimModeler4 does not have the ABAQUS 2D export):
| File > Import discrete Data > uncheck all.
| |topography in SimModeler| Now we can specify the size we want in the
  central and the side areas. Here is the resulting mesh:
| |meshed topography| The mesh can then be exported:
| File > Export mesh > ABAQUS 2D > test.inp (for example).
| We finally convert the inp file to a ts file readable by gocad using:
| ``python convertinp2ts.py test.inp --isolate``

.. |topography in SimModeler| image:: https://www.geophysik.uni-muenchen.de/~ulrich/fine2coarse2.png
.. |meshed topography| image:: https://www.geophysik.uni-muenchen.de/~ulrich/fine2coarse.png

