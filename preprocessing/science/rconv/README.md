# Rconv
Rconv is a small tool which transforms file given in the Standard Rupture Format (SRF) to the intermediate NetCDF Rupture Format (NRF) which is required by SeisSol for simulating kinematic rupture models.

## Building rconv
You need to have the proj.4 and the NetCDF libraries installed.
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
