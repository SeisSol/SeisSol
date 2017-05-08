# Rconv

Rconv is a small tool which transforms file given in the Standard Rupture Format (SRF) to the intermediate NetCDF Rupture Format (NRF) which is required by SeisSol for simulating kinematic rupture models.

## Building rconv
You need to have the proj.4 and the NetCDF libraries installed and make sure that the system is able to find them. Then just enter
`scons` in the main folder in order to compile rconv.

Here is a way to configure your installation before running scons:   
`export CC='icc'`   
`export CXX='icpc'`   
`export LIBRARY_PATH=Path-to-proj4/lib:Path-to-netcdf/lib:$LIBRARY_PATH`   
`export CPATH=Path-to-proj4/include:Path-to-netcdf/include:Path-to-SeisSol/SeisSol/src:Path-to-SeisSol/SeisSol/submodules$CPATH`

## Using rconv
Starting rconv without arguments gives you a short introduction for using the tool. You may furthermore consult the [Wiki entry](https://github.com/SeisSol/SeisSol/wiki/Standard-Rupture-Format) about the Standard Rupture Format.
