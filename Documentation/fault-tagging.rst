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

Note that ``<boundaryCondition tag="3">`` is equivalent to ``<dynamicRupture>``. Therefore, if you want to tag face 2 as 3, you can use either: 

.. code-block:: xml

   <boundaryCondition tag="3">2</boundaryCondition> 

or

.. code-block:: xml

   <dynamicRupture>2</dynamicRupture> 

Note also that if a face is tagged twice, only the first tag will be considered. 


Using more than 189 dynamic rupture tags
----------------------------------------

To handle more than 189 dynamic rupture tags (i.e. more than 255 fault tags), you will need to adjust the boundary format to at least i64 (or i32x4).

That is, build your mesh in PUMgen with the option ``--boundarytype=int64``.
Next, specify ``pumlboundaryformat = 'i64'`` in the ``&meshnml`` section of your SeisSol parameter file.

The i64 boundary format has an upper limit of 65469 dynamic rupture tags (65535 fault tags).
For an even larger tag space, you will need to resort to the i32x4 format. For that, proceed
as before, replacing ``i64`` by ``i32x4``.

To see which boundary format you have built your mesh for, you can use ``h5dump -H <yourmeshfile>.puml.h5``,
and look at the datatype and the shape of the ``boundary`` dataset.

In PUMgen, meshes with more than 255 tags can be created using ``pumgen -xml`` option, e.g. :

.. code-block:: xml

   <boundaryCondition tag="3">13245</boundaryCondition>
   .
   .
   .
   <boundaryCondition tag="900">12345,14325</boundaryCondition>


