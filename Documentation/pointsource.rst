Point Source
============

SISMOWINE WP2\_LOH1
~~~~~~~~~~~~~~~~~~~

SISMOWINE is intended as a long-term interactive web interface for
verifying numerical modeling methods in seismology. Numerical-method
developers and numerical modelers may compare their solutions with other
solutions. SISMOWINE is a continuation of the original SPICE Code
Validation interface established within the 6th Framework Programme
project .

LOH1 is used as an example here to illustrate the implementation of
source point for earthquake nucleation in SeisSol. The details of LOH1
model can also be found at .

The model uses Right-handed Cartesian, x positive North, y positive
East, z positive downward, all coordinates in meters. The source is
buried at 2000 m in a half-space Earth (Figure [fig:loh1]. The top layer
is 1000 m thick and the bottom layer is 33000 m. The material parameters
are listed in Table [table:loh1].

+--------------+------------+-----------+-----------+-------+-------+
|              | Vp (m/s)   | Vs(m/s)   | density   | Qp    | Qs    |
+==============+============+===========+===========+=======+=======+
| layer        | 4000       | 2000      | 2600      | Inf   | Inf   |
+--------------+------------+-----------+-----------+-------+-------+
| half-space   | 6000       | 3464      | 2700      | Inf   | Inf   |
+--------------+------------+-----------+-----------+-------+-------+

Table: Material properties in LOH1 .

.. figure:: LatexFigures/LOH1.jpg
   :alt: Geometry of LOH1 .
   :width: 11.00000cm

   Geometry of LOH1 .

Geometry
~~~~~~~~

The mesh is generate using Gmsh.

.. figure:: LatexFigures/loh1_mesh.png
   :alt: Geometry of LOH1 model (Gmsh)
   :width: 11.00000cm

   Geometry of LOH1 model (Gmsh). A 1 km layer of low velocity (Vp=4000
   m/s, vs=2000 m/s) is at the top of high velocity (vp=6000 m/s,
   vs=3464 m/s).

Point source input
~~~~~~~~~~~~~~~~~~

The point source needs to be turned on in *parameter.par* file.
:: 

  &SourceType
  Type = 50
  FileName=’LOH1_source.dat’
  /

The source input file can be found at . Duration of the source is 4
seconds.

Results
~~~~~~~

| The comparison with solution is shown in Figure [fig:compare\_loh1].

.. figure:: LatexFigures/loh1_benchmark.png
   :alt: Benchmark of x-component particle velocity
   :width: 11.00000cm

   Benchmark of x-component particle velocity at receiver point 1 (0.0,
   693.0,0.1). Bule is 4-order SeisSol and orange is SISMOWINE result. 

   
   
