export KMP_AFFINITY=proclist=[0-$((OMP_NUM_THREADS-1))],granularity=thread,explicit

