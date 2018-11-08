easi
====

Introduction
------------

easi is a library for the Easy Initialization of model parameters in
three (or less) dimensional domains. easi offers the possibility to
parameterize the simulation without having to recompile SeisSol. Thanks
to easi, all user can run their simulations with the same executable (no
more hardcoded fortran initialization routines).

Easi Documentation
------------------

| A general description of easi can be found on the README:
| `https://github.com/SeisSol/easi <https://github.com/SeisSol/easi>`__

| A detailed description of all easi component can be found in the wiki:
| `https://github.com/SeisSol/easi/wiki <https://github.com/SeisSol/easi/wiki>`__

Debugging easi script
---------------------

| Most easi components return easy to track error, for example
| ``test.yaml: yaml-cpp: error at line 6, column 9: illegal map value``
| Yet implajit function map are more complex to debug. The following
  example:
| ``27.1: syntax error, unexpected '}', expecting ;``
| indicates that an error occur in the 27th line of the function, but
  does not indicate which file and which function.
| Hopefully this will be improved in the future.
