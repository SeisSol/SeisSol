## Source Code Documentation
This documentation is intended to be an extra tool for people involved in SeisSol development to help to navigate within the source code. There are two types of the documentation, namely: a short one, for a quick reference, and an extended one, which includes some information from SeisSol submodules. 

### Getting Started
Please, make sure that you have Doxygen installed
```concole
$ apt-get install doxygen
```

Don't forget to download SeisSol submodules in case if you want to generate the extended documentation.

```console
$ cd ..
$ git submodule update --init --recursive
$ cd SourceCodeDocs 
```

### Generation

Execute the following command in order to generate and view the short form:
```
$ doxygen Doxyfile
$ firefox ./docs/html/index.html
```
or the extended one:
```
$ doxygen ExtendedDoxyfile
$ firefox ./extended_docs/html/index.html
```
##### Note
You can use your favorite browser instead of *firefox* to open a documentation
