Initial Conditions
==================

Currently we provide the following initial conditions:

Zero        
----

All quantities are set to zero. 
This is the standard case to work with point sources or dynamic rupture.

Planarwave  
----------

A planar wave for convergence tests.
The inital values are computed such that a planar wave in a unit cube is imposed.
This scenario needs periodic boundary conditions to make sense.
This is the only case where the old netcdf mesh format is prefered.
After the simulation is finished the errors between the analytic solution and the numerical one are plotted in the :math:`L^1`-,  :math:`L^2`- and :math:`L^\infty`-norm.

Use ``cube_c`` to generate the meshes for the convergence tests:
https://github.com/SeisSol/SeisSol/tree/master/preprocessing/meshing/cube_c

Superimposed Planarwave
-----------------------

Superimposed three planar waves travelling in different directions. 
This is especially interesting in the case of directional dependent properties such as for anisotropic materials.

Scholte     
-------

A Scholte wave to test elastic-acoustic coupling

Snell       
-----

Snells law to test elastic-acoustic coupling

Ocean       
-----

An uncoupled ocean test case for acoustic equations


How to implement a new initial condition?
-----------------------------------------

New initial conditions can be easily implemented. Extend the class 

.. code-block:: c

  seissol::physics::Initalfield

and implement the method 

.. code-block:: c

  void evaluate(  double time,
                  std::vector<std::array<double, 3>> const& points,
                  const CellMaterialData& materialData,
                  yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const;

Here :code:`dofsQP(i,j)` is the value of the :math:`j^\text{th}` quantity at the :code:`points[i]`.
