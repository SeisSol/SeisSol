  
# Description

This folder contains scripts to compute Ground motions parameters
(e.g. PGA, PGV, PGD, SA(T)) from a SeisSol surface output file.

# Prerequisites

The main script (ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py)
requires gmpe-smtk. The latest gmpe-smtk depends on openquake, but some older
versions do not require this dependency. Therefore, one can choose to either
install the deprecated gmpe-smtk without openquake, or the latest gmpe-smtk with
openquake.

## Installing a deprecated gmpe-smtk

In this case, no installation is needed, the repository only needs to be cloned:

```bash
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

<!-- markdownlint-disable MD013 -->
```bash
pip install -r https://github.com/gem/oq-engine/raw/master/requirements-py311-linux64.txt openquake.engine
pip install openquake.engine
```
<!-- markdownlint-enable MD013 -->

On SupermucNG, see <https://seissol.readthedocs.io/en/latest/supermuc.html#accessing-pypi>
for using pip.

The installation on other operating system that Linux is documented here:
<https://github.com/gem/oq-engine/tree/master/doc/installing>

You can test the installation by running an example:

```bash
oq engine --help
#alternative
git clone https://github.com/gem/oq-engine.git --depth 1
oq run oq-engine/demos/hazard/AreaSourceClassicalPSHA/job.ini
```

Then the latest gmpe-smtk can be cloned and added to the pythonpath as previously
described.

## Installing other python modules

Simply use:

```bash
pip install -r requirements.txt
```

If not root, use ``--user``.

# Running the script

All the available option for running the script can be displayed using:

```bash
python3 ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py -h
```

A typical example for running the script on a desktop can be:

<!-- markdownlint-disable MD013 -->
```bash
python3 ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py --MP 4 prefix-surface.xdmf --noMPI
```
<!-- markdownlint-enable MD013 -->

(assuming a 4-core processor).

The script can also be used on a high performance cluster, see
`SupermucNGCommandFileComputeGroundMotions.sh` for an example of a job file on supermucNG.
