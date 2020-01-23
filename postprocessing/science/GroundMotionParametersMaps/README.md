
# Description

This folder contains scripts to compute Ground motions parameters (e.g. PGA, PGV, PGD, SA(T)) from a SeisSol surface output file.

# Prerequisites

The main script (ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py) requires gmpe-smtk. No installation is needed, the repository only needs to be cloned:

```
git clone https://github.com/GEMScienceTools/gmpe-smtk
cd gmpe-smtk
git checkout 4f008173e89f6e4ba4450fb43e95ffdf51a7c2ba^
#now add gmpe-smtk to your python path
export PYTHONPATH=$PYTHONPATH:'$(pwd)
##you could also consider adding it to your ~.bashrc
#echo 'export PYTHONPATH=$PYTHONPATH:'$(pwd) >> ~/.bashrc

```

In addition, you may need to install missing python modules, e.g. for mpi4py (if running the script on multiple ranks). This can be done e.g. through anaconda, using:

```
conda install -c anaconda mpi4py
```
On a machine without internet, this can be done through pip. First download the module on a computer with access to internet (example for module lxml):

```
pip download lxml
```

then transfer the \*.whl to the cluster. and install the module using:
```
pip install --user lxml-4.4.1-cp36-cp36m-manylinux1_x86_64.whl
```



# Running the script
 
All the available option for running the script can be displayed using:

```
python ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py -h
```

A typical example for running the script on a desktop can be:
```
python ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py --MP 4 prefix-surface.xdmf --noMPI
```
(assuming a 4 cores processors).

The script can also be used on a high performance cluster, see Supermuc2CommandFileComputeGroundMotions.sh and SupermucNGCommandFileComputeGroundMotions.sh for examples of Job File on supermuc2 and supermugNG.
