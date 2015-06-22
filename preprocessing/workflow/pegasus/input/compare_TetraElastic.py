#!/bin/env python
##
# @file
# This file is part of SeisSol.
#
# @author Fabio Gratl (f.gratl AT in.tum.de)
#
# @section LICENSE
# Copyright (c) 2013, SeisSol Group
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# @section DESCRIPTION
# Takes 2 stdout files
# Extracts MPI Error analysis
# compares difference
# print everything

import sys

#load reference File
r_name=sys.argv[1]
r=open(r_name)
#List for reference values
refList=[]
#search reference file for values
for line in r:
        if line.find("MPI Error analysis") > -1:
                rvar = line.split()[-1]
                if not(rvar in  refList):
				  for line in r:
						  if line.find("L1_norm") > -1:
								  rvalL1 = line.split()[-1]
								  #print val
						  if line.find("L2_norm") > -1:
								  rvalL2 = line.split()[-1]
						  if line.find("Linf_norm") > -1:
								  rvalLinf = line.split()[-1]
								  rtouple = rvar, rvalL1, rvalL2, rvalLinf
								  #print touple
								  refList.append(rtouple)
								  break
print ""
print "Reference Values from File:"
print r_name
print "variable \t\tL1 norm \t\t\tL2 norm \t\t\tL inf norm"
for line in refList:
  for tRefList in line:
	print "%s             \t" %(tRefList),
  print

#load output File
f_name=sys.argv[2]
f=open(f_name)

#search output File for values
valueList=[]
for line in f:
	if line.find("MPI Error analysis") > -1:
		var = line.split()[-1]
		#print var
		for line in f:
			if line.find("L1_norm") > -1:
				valL1 = line.split()[-1]
				#print val
			if line.find("L2_norm") > -1:
				valL2 = line.split()[-1]
			if line.find("Linf_norm") > -1:
				valLifnf = line.split()[-1]
				touple = var, valL1, valL2, valLifnf
				#print touple
				valueList.append(touple)
				break
		
print ""
print "Benchmark Values from File:"
print f_name
print "variable \t\tL1 norm \t\t\tL2 norm \t\t\tL inf norm"
for line in valueList:
  for tValueList in line:
	print "%s             \t" %(tValueList),
  print 

#List with differences
diffList=[]
for refTouple in refList:
	for valTouple in valueList:
		if refTouple[0] == valTouple[0]:
			diffL1 = float(valTouple[1]) - float(refTouple[1])
			diffL2 = float(valTouple[2]) - float(refTouple[2])
			diffLinf = float(valTouple[3]) - float(refTouple[3])
			diffTouple = refTouple[0], diffL1, diffL2, diffLinf
			diffList.append(diffTouple)
	
print ""
print "Difference"
print "variable \t\tL1 norm \t\t\tL2 norm \t\t\tL inf norm"
for line in diffList:
	for t in line:
	  print "%s             \t" %(t),
	print
