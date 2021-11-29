IDEs for development
====================

Clion
-----

We highly recommend to use `Clion` for your development. It has a robust support for C/C++ languages and `CMake`.


Eclipse for Fortran
-------------------

Unfortunately, Clion does not support Fortran well. Therefore, we recommend to have an installation for eclipse parallel edition which comes with a Photran plugin.

Build SeisSol as usual with CMake but from the project root directory (no from build folder). You will need to extend your cmake commands during a project configuration

```
cmake . <settings> -G"Eclipse CDT4 - Unix Makefiles"
make
```

Once it is done, you will need to an a `fortran nature` to the project description. Edit the nature block in `.project` file as follow:

```
	<natures>
		<nature>org.eclipse.cdt.make.core.makeNature</nature>
		<nature>org.eclipse.cdt.make.core.ScannerConfigNature</nature>
		<nature>org.eclipse.cdt.core.ccnature</nature>
		<nature>org.eclipse.cdt.core.cnature</nature>
		<nature>org.eclipse.photran.core.fnature</nature>
	</natures>

```

Open Eclipse and the project. Go to `Project/Properties/Fortran General/Analysis and Refactoring` and click on `Enable Fortran analysis/refactoring`, and click `Apply`. You will need to close and open Eclipse, and build the project.

Note, that Photran can detect some minor syntactical errors which are not crucial for GNU and Intel compilers i.e., new line as `\` instead of `&`. Fix them! It is displayed on the right panel as `Syntax Exceptions`. Otherwise, searching declarations is not going to work. Once it is done, re-compiler the code.

It will take some time for Photran to analyze the Fortran source code and build its internal data structures. During this time it seems that Eclipse doing nothing. Wait, it will come!