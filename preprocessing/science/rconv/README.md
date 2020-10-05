# rconv
rconv is a tool which allows converting files describing kinematic rupture models from the Standard Rupture Format (*.srf) to the intermediate NetCDF Rupture Format (*.nrf).

## Building rconv

rconv depends on the PROJ.4 and NetCDF libraries.

### Installing PROJ.4

rconv relies on a deprecated version of the PROJ library (as it includes the file 'projects.h', which is not anymore built in recent releases of PROJ, see [OSGeo/PROJ#835](OSGeo/PROJ#835)). Also, the most recent version of PROJ need access to the internet when running cmake (because of the googletest 1.8.0 framework). To install PROJ on supermuc (behind a firewall), we therefore checkout an deprecated version of the PROJ library:

```
module load gcc
git clone git@github.com:OSGeo/PROJ
cd PROJ
git checkout 4.9.3
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$my_proj4_install_prefix
make 
make install
```
Note that `$my_proj4_install_prefix` should be different than the build directory (else `make install` will raise the error `cannot find libproj.so.12.0.0`).

### Installing netcdf

see this [link](https://seissol.readthedocs.io/en/latest/compilation.html#installing-netcdf)

### Building rconv

To install rconv with cmake, execute the following in the main folder:

```
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=$my_proj4_install_prefix
make 
```

Note that rconv can also be built with scons, using:

```
scons compiler=your-compiler netcdfDir=$NETCDF_DIR proj4Dir=Path-to-proj4
```

## Using rconv

Make sure to add this line to ~/.bashrc file:
`export LD_LIBRARY_PATH=path-to-proj.4/build/lib:$LD_LIBRARY_PATH`
Starting rconv without arguments gives you a short introduction for using the tool. You may furthermore consult the [Documentation](https://seissol.readthedocs.io/en/latest/standard-rupture-format.html) about the Standard Rupture Format.

## Dealing with projected data

If the SRF data are already projected, the projection within rconv can be by-passed by compiling rconv with:
`scons compiler=your-compiler netcdfDir=$NETCDF_DIR proj4Dir=Path-to-proj4 NoProj=True`
