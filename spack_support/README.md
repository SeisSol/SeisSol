# Spack and SeisSol
Installing any HPC software can be tricky. Unfortunately, SeisSol is not an exception. To considerably alleviate the installation process, we provide few scripts which can be integrated into [Spack](https://github.com/spack/spack/wiki) which is a new HPC software package manager. The installation processes is devided into two parts, namely: **seissol-env** and **seissol-core**.

**seissol-env**  installs all libraries that SeisSol depends on for a given compiler. For instance, mpi, hdf5, netcdf, asagi, etc. The scrip targets developers and users who are used to, or have to, intsall SeisSol manually.

**seissol-core** installs SeisSol itself. It depends and invokes **seissol-env**, however, it is intended for the end users who don't want to bother themselves with manual compilation of SeisSol with either scons and cmake.


## Getting Started
First of all, you have to install Spack which you can download from the official github [repo](https://github.com/spack/spack.git). Make sure that you have [git](https://en.wikipedia.org/wiki/Git) and python installed on your system.

```console
:~$ cd $HOME
:~$ git clone https://github.com/spack/spack.git
```
Append you *.bashrc* to have Spack avaliable all the time
```console
:~$ cd spack
:~$ echo "export SPACK_ROOT=$PWD" >> $HOME/.bashrc
:~$ echo "export PATH=\$SPACK_ROOT/bin:\$PATH" >> $HOME/.bashrc
:~$ cd $HOME
```
Close and open your terminal to make sure that the changes have been applied. For the next step, you need to acquire Spack scripts for SeisSol. Download SeisSol from [here](https://github.com/SeisSol/SeisSol) and go to the root directory of the application
```console
:~$ git clone https://github.com/SeisSol/SeisSol.git
:~$ cd SeisSol
```
To make SeisSol installation scripts be visiable inside Spack, one has to copy them to Spack as following:
```console
:~$ cp -rf spack_support/seissol-* $SPACK_ROOT/var/spack/repos/builtin/packages/
```
To make sure that everything went well, query avaliable packages in Spack.
```console
:~$ spack list seissol*
==> 2 packages.
seissol-core  seissol-env
```
If you can see the output similar as above then we are ready to proceed!

### Prerequisites
One of the main ideas of Spack is to produce a consistent build of your software stack, i.e. when everything is compiled with the same set of compilers. You may have your preferable compilers installed on your system. If so, you can add them to Spack.
```console
:~$ spack compiler find <path_to_your_compiler>
```
However, if you don't have any or you want to try another one you can install it with Spack. 
For example, you have to do the following to install gcc 8.3.0:
```console
spack install gcc%8.3.0
```
Don't forget to add it to Spack once it has been installed:
```console
:~$ spack compiler find $(spack location -i gcc@8.3.0)
```


### Install options
##### seissol-env
As it was mentioned above, this script is responsible for installing all libs that SeisSol depends on. By default, the script tries to install the following packages: *netcdf-c-4.6.1*, *hdf5-1.8.21*, *libxsmm*, *metis (with 64bit-wide integer support)*, *memkind*. Additionally, some options, that we are going to discuss below, can triger installing additional libraries. For example, *parmetis*, *intel mkl*, *asagi*, *mpi*, etc.

The script has the following options:
| Options       | Description                              |
| ------------- |:----------------------------------------:|
| mpi      | configures all libs to provide mpi support    |
| openmp   | configures all libs to provide multithreading |
| asagi    | installs Asagi                                |
| blas     | installs a blas implementation                |

Note, that some flags are in charge of so-called "virtual" packages like *mpi* and *blas*, and they have to be clarified. 

Let's assume that we have Intel 17.0.2 compiler and we want to isntall all default libs and Asagi with MPI support. Additionally, we want to build (if it hasn't been built before) and compile eveything with OpenMPI (version 3.1.5) implemetation of the MPI standard. Then we have to execute the following:
```console
:~$ spack install seissol-env +mpi +asagi %intel@17.0.2 ^openmpi@3.1.5
```
As you can see, we have to clarify what MPI implementation we want to use. The same for blas.
```console
:~$ spack install seissol-env +mpi +blas +asagi %intel@17.0.2 ^openmpi@3.1.5 ^openblas@0.3.4
```
Please, refer to the Spack [documentation](https://spack.readthedocs.io/en/latest/) to learn more about how to influence and provide build options for dependencies. For instance, if one wants to install and configure SeisSol software stack with Cuda Aware MPI, the following should be executed:
```console
:~$ spack install seissol-env +mpi %gcc@8.3.0 ^openmpi@3.1.5+cuda ^metis+int64 ^cuda@10.1.243
```
Additionally, one can notice what we require either to use or to install *CUDA 10.1*

If you have some packages installed on your system and you don't want Spack to install versions of these packages again and again, you need to edit *packages.yaml* file
```console
:~$ vi ~/.spack/packages.yaml
```
Please, refere to the Spack [documentation](https://spack.readthedocs.io/en/latest/build_settings.html#build-settings) to learn how to configure the file.

##### seissol-core
Blank. It is not ready
### Workflow
blank
### Examples
blank
### Known isssues
blank
