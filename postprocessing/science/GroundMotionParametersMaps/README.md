  
# Description

This folder contains scripts to compute Ground motions parameters (e.g. PGA, PGV, PGD, SA(T)) from a SeisSol surface output file.

# Prerequisites

The main script (ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py) requires gmpe-smtk. The latest gmpe-smtk depends on openquake, but some older versions do not require this dependency. Therefore, one can choose to either install the deprecated gmpe-smtk without openquake, or the latest gmpe-smtk with openquake.

## Installing a deprecated gmpe-smtk
In this case, no installation is needed, the repository only needs to be cloned:

```
git clone https://github.com/GEMScienceTools/gmpe-smtk
cd gmpe-smtk
git checkout 4f008173e89f6e4ba4450fb43e95ffdf51a7c2ba^
#now add gmpe-smtk to your python path
export PYTHONPATH=$PYTHONPATH:'$(pwd)
##you could also consider adding it to your ~.bashrc
#echo 'export PYTHONPATH=$PYTHONPATH:'$(pwd) >> ~/.bashrc

```
## Installing the latest gmpe-smtk and openquake

Here is a proposed workflow for installing openquake:
```
#Clone the repository
git clone https://github.com/gem/oq-engine.git
#Use pip to install and specify dir_to_install with --target
#This step requires python 3.6 and some other packages
pip install -r oq-engine/requirements-py36-linux64.txt -r oq-engine/requirements-extra-py36-linux64.txt --target dir_to_install
#add the oq-engine to your pythonpath
export PYTHONPATH=$PYTHONPATH:/dir_to_qo-engine/
```
The installation on other operating system that Linux is documented here:
https://github.com/gem/oq-engine/blob/master/doc/installing/development.md

You can test the installation by running an example:

```
oq-engine --run dir_to_oq-engine/demos/hazard/AreaSourceClassicalPSHA/job.ini
#alternative
oq engine --help
```
Then the latest gmpe-smtk can be cloned and added to the pythonpath as previously described.

## Installing other python modules

Depending on your system, you may need to install missing python modules, e.g. for mpi4py (if running the script on multiple ranks). This can be done e.g. through anaconda, using:

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

The script can also be used on a high performance cluster, see SupermucNGCommandFileComputeGroundMotions.sh for an example of a job file on supermugNG.
