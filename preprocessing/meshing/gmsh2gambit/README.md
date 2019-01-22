gmsh2gambit
===============
Transform msh file (gmsh msh 2.0 mesh file) to neu file (neutral mesh file)
From gmsh 4.0, the option '-format neu' of gmsh should be used instead.

Compilation
===============
```
export CC=gcc (or icpc with intel)
export CXX = g++ (or icc with intel)
scons   
```

Use
===============
gmsh2gambit -i test.msh -o test.neu


