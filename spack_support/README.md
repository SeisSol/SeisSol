# SeisSol and Spack v. 0.13.3
Installing any HPC software can be tricky. Unfortunately, SeisSol is not an exception. To considerably alleviate the installation process, we provide few scripts which can be integrated into [Spack](https://github.com/spack/spack/wiki) which is a new HPC software package manager. The installation processes is devided into two parts, namely: **seissol-env** and **seissol-core**.

**seissol-env**  installs all libraries that SeisSol depends on for a given compiler. For instance, mpi, hdf5, netcdf, asagi, etc. The scrip targets developers and users who are used to, or have to, intsall SeisSol manually.

**seissol-core** installs SeisSol itself. It depends and invokes **seissol-env**, however, it is intended for the end users who don't want to bother themselves with manual compilation of SeisSol with either scons and cmake.


## Getting Started
First of all, you have to install Spack which you can download from the official github [repo](https://github.com/spack/spack.git). Make sure that you have [git](https://en.wikipedia.org/wiki/Git), [curl](https://curl.haxx.se/libcurl/) and python installed on your system. Also, be sure that your Linux distribution has *build-essential* packages.

You can install it as following:
```console
:~$ apt-get install git python curl build-essential
```
However, most of Unix-kind operating systems come with these packages pre-installed. Now, let's install Spack.
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
To make SeisSol installation scripts be visiable inside of Spack, one has to add them to the Spack repository. We recomend to install our scripts into a separete directory to avoid problems with dangling files inside of Spack in case if you decide to delete the current SeisSol repository.

```console
$ cd spack_support
$ mkdir build && cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=<install_dir>
$ spack repo add <install_dir>/spack_support
```
To make sure that everything went well, query avaliable packages in Spack.
```console
:~$ spack list seissol*
==> 2 packages.
seissol-core  seissol-env
```
If you can see the output similar as above then we are ready to proceed!

Please, keep in mind that we update installation scripts from time to time. Therefore, you have remove old scripts from spack as following:
```console
:~$ spack repo remove spack_support
```

Don't forget to add new ones into the spack in the same way how we did above.

### Prerequisites
One of the main ideas of Spack is to produce a consistent build of your software stack, i. e. when everything is compiled with the same set of compilers. You may have your preferable compilers installed on your system. If so, you can add them to Spack.
```console
:~$ spack compiler find <path_to_your_compiler>
```
However, if you don't have any or you want to try another one you can install it with Spack. 
For example, you have to do the following to install gcc 8.3.0:
```console
spack install gcc@8.3.0
```
Don't forget to add it to Spack once it has been installed:
```console
:~$ spack compiler find $(spack location -i gcc@8.3.0)
```

### Install options
##### seissol-env
As it was mentioned above, this script is responsible for installing all libs that SeisSol depends on. By default, the script tries to install the following packages: *netcdf-c-4.6.1*, *hdf5-1.8.21*, *libxsmm*, *metis (with 64bit-wide integer support)*, *memkind*. Additionally, some options, that we are going to discuss below, can triger installing additional libraries. For example, *parmetis*, *intel mkl*, *asagi*, *mpi*, etc.

The script has the following options:

| Options  | Descriptions                                  |
| ---------|:---------------------------------------------:|
| mpi      | configures all libs to provide mpi support    |
| openmp   | configures all libs to provide multithreading |
| asagi    | installs Asagi                                |
| blas     | installs a blas implementation                |

Note, that some flags are in charge of so-called "virtual" packages like *mpi* and *blas*, and they have to be clarified. 

Let's assume that we have Intel 17.0.2 compiler and we want to install all default libs and Asagi with MPI support. Additionally, we want to build (if it hasn't been built before) and compile everything with OpenMPI (version 3.1.5) implementation of the MPI standard. Then we have to execute the following:
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
###### How can I use these installed packages?
Fortunately, Spack allows you to install environment-modules to manage your enviroment variables
```console
:~$ spack bootstrap
:~$ source $SPACK_ROOT/share/spack/setup-env.sh
```
You can add the last command to your *bashrc* script to be able to access installed modules whenever you start a new terminal. However, it might result in increasing start-up time of your terminal.
```console
:~$ echo "source \$SPACK_ROOT/share/spack/setup-env.sh" >> $HOME/.bashrc
```
After that you can access installed libraries using *module avail*, *module list*, *module load*, etc. 
##### seissol-core
Blank. It is not ready

### Examples
Seissol with GNU tools:
```console
:~$ spack install seissol-env +mpi +asagi %gcc@8.3.0 ^openmpi@3.1.5 ^metis+int64
```

Seissol with Intel tools:
```console
:~$ spack install seissol-env +mpi +asagi %intel@17.0.2 ^intel-mpi@2018.0.128 ^metis+int64
```
### Known issues
1. Spack is quite sensitive with respect to your environment. It can mistakenly fetch/reference a package that you installed manually and added to one of the following env. variables: PATH, LD_LIBRARY_PATH, LIBRARY_PATH, CPATH, CMAKE_PREFIX_PATH, PJG_CONFIG_PATH, etc. Please, try to keep these variables as "clean" as possible and leave only essentials.

2. Even though you removed all old manually installed libraries/packages from your env. variables to keep your environment "clean", be sure that the env. variables do not contain syntactical mistakes. The example below shows a small syntactical mistake in *$LD_LIBRARY_PATH* but it does not usually cause any problem for most of operating systems. However, the dot in the middle can break compilation of some packages under Spack.
```console
:~$ echo $LD_LIBRARY_PATH
<path_1>:.:<path_2>
```


3. Some compilers, especially new ones, are not always able to successfully install all SeisSol software stack. The solution is to look at the building log-file of Spack. If you suspect a compiler issue then try to use a previous version of your compiler and try again.

4. ImpalaJIT has a limitted support w.r.t GNU compilers. It compiles with gcc v. 5.5.0 or lower.