# rconv

rconv is a tool which allows converting files describing kinematic rupture models
from the Standard Rupture Format (\*.srf) to the intermediate
NetCDF Rupture Format (\*.nrf).

## Installing rconv

rconv depends on the PROJ.4 and NetCDF libraries.

### Installing PROJ.4

rconv relies on a deprecated version of the PROJ library (as it includes the file
`projects.h` which is not anymore built in recent releases of PROJ, see
<https://github.com/OSGeo/PROJ/issues/835>). Also, the most recent version of
PROJ needs access to the internet when running cmake
(because of the googletest 1.8.0 framework).
To install PROJ on supermuc (behind a firewall), we therefore
checkout a deprecated version of the PROJ library:

```bash
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

Note that `$my_proj4_install_prefix` should be different than the build directory
(else `make install` will raise the error `cannot find libproj.so.12.0.0`).

### Installing netcdf

see this [link](https://seissol.readthedocs.io/en/latest/compilation.html#installing-netcdf)

### Building rconv

To install rconv, execute the following in the main folder:

```bash
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=$my_proj4_install_prefix
make 
make install
```

## Using rconv

Starting rconv without arguments gives you a short introduction for using the tool.
You may furthermore consult the [Documentation](https://seissol.readthedocs.io/en/latest/standard-rupture-format.html)
about the Standard Rupture Format.

## Compilation without PROJ.4

rconv can be compiled without proj (e.g. by dropping the argument `-DCMAKE_PREFIX_PATH=$my_proj4_install_prefix`)
In that case, the input file won't be projected anymore.
