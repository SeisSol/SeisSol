.. _Configuration:

Configuration
=============

To set up a SeisSol run, you need to do the following steps (assuming,
your are in the root directory of the code repository):

1. Create the launch directory:

   ::

       mkdir launch_SeisSol

2. Copy all executables to the launch directory:

   ::

       cp build/SeisSol* launch_SeisSol/

3. Write path to the Maple directory to a file called ``DGPATH``:

   ::

       echo $PWD/Maple/ > launch_SeisSol/DGPATH

4. Create your :doc:`parameter-file`

5. Copy any additional input files referenced in the parameter file (for
   example file with receiver coordinates) to your launch directory

6. (For large mesh only) Create symbolic links (ln -s) to the mesh
   file(s) in your launch directory

7. Make sure output and checkpoint directories exist

8. Optional: set :doc:`environment-variables` for tuning

Checklist for required files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Necessary files
^^^^^^^^^^^^^^^

-  SeisSol executable (compiled on the system where the job will run)
-  DGPATH
-  Parameter file
-  \*.yaml files for setting model parameters

Optional files depending on settings in the parameter file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  receiver files in \*.dat format (if nRecordPoints >0 in the parameter
   file)
-  fault receiver files in \*.dat format (in the parameter file)
