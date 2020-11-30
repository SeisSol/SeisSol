Fault tagging
=============

| In SeisSol, boundary conditions are tagged as:
| 0: regular
| 1: free surface
| 2: free surface + gravity (water surface)
| 3 or n>64: dynamic rupture
| 5: absorbing
| 6: periodic

Dynamic rupture can therefore be tagged using a range of possible tags.
This allows initializing fault parameters segment-wise
easily. For example, if we have 2 segments, and we want them to have
different dynamic friction, we can tag them with 3 and 65 and then use:

.. code-block:: yaml

   [mu_d]: !Any
     components:
       - !GroupFilter
         groups: 3
         components: !ConstantMap
           map:
             mu_d:    0.3
       - !GroupFilter
         groups: 65
         components: !ConstantMap
           map:
             mu_d:    0.4

Currently, the only way to tag fault faces other tags than 3 with SimModeler is to use the `--xml` option of pumgen. 
For example, to tag face 2 as 3 and face 8 and 9 as 65, we would
use:

.. code-block:: xml

   <boundaryCondition tag="3">2</boundaryCondition>
   <boundaryCondition tag="65">8,9</boundaryCondition>

Then pumgen is run using the xml option:

::

   pumgen -s simmodsuite -l SimModelerLib.lic --xml MeshandAnalysisAttributes.xml prefix.smd output_prefix

