# Spack and SeisSol
Installing any HPC software can be a tricky task. Unfortunately, SeisSol is not an exception. To considerably alleviate the installation process, we provide few scripts which can be integrated to [Spack](https://github.com/spack/spack/wiki) which is a new HPC software package manager. The installation processes is devided into two parts, namely: **seissol-env** and **seissol-core**.

**seissol-env**  installs all libraries that SeisSol depends on for a given compiler. For instance, mpi, hdf5, netcdf, asagi, etc. The scrip targets developers and users who are used to, or have to, intsall SeisSol manually.
**seissol-core** installs SeisSol itself. It depends and invokes **seissol-env**, however, it is intended for the end users who don't want to bother themselves with manual compilation of SeisSol with either scons and cmake.


## Getting Started
First of all, you have to install Spack which you can download from the official github [repo](https://github.com/spack/spack.git). 

### Installation options
blank
### Workflow
blank
### Examples
blank
### Known isssues
blank
