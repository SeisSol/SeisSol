Fault tagging
=============

| In SeisSol, boundary conditions are historically tagged as:
| 0: regular
| 1: free surface
| 3: dynamic rupture
| 5: absorbing
| 6: periodic

Recently, we added the possibility to also tag dynamic rupture also with
n>6. It is then possible to initialized fault parameters segment-wise
easily. For example, if we have 2 segments, and we want them to have
different dynamic friction, we can tag them with 3 and 10 and then use:

.. code-block:: yaml

   [mu_d]: !Any
     components:
       - !GroupFilter
         groups: 3
         components: !ConstantMap
           map:
             mu_d:    0.3
       - !GroupFilter
         groups: 10
         components: !ConstantMap
           map:
             mu_d:    0.4

Currently the only way to tag fault faces using simModeler library and
pumgen is to compile pumgen in the xml branch and make use of the xml
feature. For example to tag face 2 as 3 and face 8 and 9 as 15, we would
use:

.. code-block:: xml

   <boundaryCondition tag="3">2</boundaryCondition>
   <boundaryCondition tag="15">8,9</boundaryCondition>

Then pumgen is run using the xml option:

::

   poe pumgen -s simmodsuite -l SimModelerLib.lic --xml MeshandAnalysisAttributes.xml prefix.smd output_prefix

