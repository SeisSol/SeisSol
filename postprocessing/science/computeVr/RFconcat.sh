#!/bin/sh
#Concatenate all the RF files into a single file
folder="Resultstpv29Z192b"
prefix="tpv29Z"

cat $folder/$prefix-RF-* > trash.dat
awk 'NF==1 {sum=sum+$1} END {print sum}' trash.dat > $prefix-RF-concat.dat
awk 'NF==4 {print $0}' trash.dat >> $prefix-RF-concat.dat
rm trash.dat
