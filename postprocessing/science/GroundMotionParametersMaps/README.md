# Description

This folder contains scripts to compute Ground motions parameters (e.g. PGA, PGV, PGD, SA(T)) from a SeisSol surface output file.

# Prerequisites

The main script (ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py) requires gmpe-smtk. No installation is needed, the repository only needs to be cloned:

```
git clone https://github.com/GEMScienceTools/gmpe-smtk
cd gmpe-smtk
git checkout 4f008173e89f6e4ba4450fb43e95ffdf51a7c2ba^
```

In addition, you may need to install missing python modules, e.g. for mpi4py (if running the script on multiple ranks). This can be done e.g. through anaconda, using:

```
conda install -c anaconda mpi4py
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

The script can also be used on a high performance cluster, see Supermuc2CommandFileComputeGroundMotions.sh for an example of Job File on supermuc2.
