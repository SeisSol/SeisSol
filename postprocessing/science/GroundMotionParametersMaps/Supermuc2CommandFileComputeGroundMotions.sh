#!/bin/bash

testcasefolder=
#@ initialdir =

# DO NOT USE environment = COPY_ALL
#@ job_type = parallel
#@ class = test

#@ island_count = 1
#@ node = 20
#@tasks_per_node = 1

#@ wall_clock_limit = 0:30:00
#@ network.MPI = sn_all,not_shared,us

#@ energy_policy_tag = di73yeq_ept
#@ minimize_time_to_solution = yes

#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification=always
#@ notify_user=yourmail
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh

export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=28
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

. /etc/profile
. /etc/profile.d/modules.sh
module unload python mpi.ibm
module load python/2.7_anaconda_mpi
poe python TuSeisSolScripts/onHdf5/ComputeGroundMotionsEstimatesFromSurfaceMPI.py SimulationsResults-NZ/RES-NZ-easi5_090318-29mio/NZ-surface.xdmf --MP 28
#if you dont have mpi4py installed, run on 1 node
#module load python
#poe python TuSeisSolScripts/onHdf5/ComputeGroundMotionsEstimatesFromSurfaceMPI.py SimulationsResults-NZ/RES-NZ-easi5_090318-29mio/NZ-surface.xdmf --MP 28 --noMPI
