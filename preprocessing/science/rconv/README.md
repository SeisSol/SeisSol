# Rconv
Rconv is a small tool which transforms file given in the Standard Rupture Format (SRF) to the intermediate NetCDF Rupture Format (NRF) which is required by SeisSol for simulating kinematic rupture models.

## Building rconv

rconv depends on the PROJ.4 and NetCDF libraries.

### Installing PROJ.4

rconv relies on a deprecated version of the PROJ library (as it includes the file 'projects.h', which is not anymore built in recent releases of PROJ, see [OSGeo/PROJ#835](OSGeo/PROJ#835)). Also, the most recent version of PROJ need access to the internet when running cmake (because of the googletest 1.8.0 framework). To install PROJ on supermuc (behind a firewall), we therefore checkout an deprecated version of the PROJ library:

```
module load gcc
git clone git@github.com:OSGeo/PROJ
git checkout 4.9.3
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$(pwd)
make 
make install
```

### Installing netcdf

see this [link](https://seissol.readthedocs.io/en/latest/compilation.html#installing-netcdf)

### Building rconv

Then enter
`scons compiler=your-compiler netcdfDir=$NETCDF_DIR proj4Dir=Path-to-proj4`
in the main folder in order to compile rconv.

## Using rconv
Make sure to add this line to ~/.bashrc file:
`export LD_LIBRARY_PATH=path-to-proj.4/build/lib:$LD_LIBRARY_PATH`
Starting rconv without arguments gives you a short introduction for using the tool. You may furthermore consult the [Documentation](https://seissol.readthedocs.io/en/latest/standard-rupture-format.html) about the Standard Rupture Format.

## Dealing with projected data
If the SRF data are already projected, the projection within rconv can be by-passed by compiling rconv with:
`scons compiler=your-compiler netcdfDir=$NETCDF_DIR proj4Dir=Path-to-proj4 NoProj=True`
