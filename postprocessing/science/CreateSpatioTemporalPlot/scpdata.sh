#!/bin/sh

myruns="RF10-25 RF25-100 RF100-600 RF-all"
myruns="RF25-100 RF100-600 RF-all"

for run in $myruns
do
echo $run
mkdir $run
scp super:/gss/scratch/h019z/di73yeq/OutputPar$run/Par$run-faultreceiver* $run
done


