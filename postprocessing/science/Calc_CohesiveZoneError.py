##
# @file
# This file is part of SeisSol.
#
# @author Stephanie Wollherr (wollherr AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/ulrich)
#
# @section LICENSE
# Copyright (c) 2005-2018, SeisSol Group
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

#Description:
#script to calculate the estimated error given a discretization, order and cohesive zone width

import numpy as np

#INPUT: PLEASE CHANGE! 
#example from Landers simulation
import argparse
parser = argparse.ArgumentParser(description='compute expected errors on a range of fault output, using Wollherr et al. (2018) convergence study')
parser.add_argument('--pla',  dest='pla', action='store_true' , default = False, help='use if Plasticity is enabled in the model')
parser.add_argument('order',  help='order of accuracy used' ,type=int)
parser.add_argument('dx',  help='on-fault mesh size (m)' ,type=float)
parser.add_argument('coh',  help='minimum cohesive zone width (in m) (estimated in paraview using (RF-DS)*Vr, or using the python script estimateMinimumCohesiveZone.py. Very important, the on-fault mesh size should be small enough to resolve it!!!)' ,type=int)
args = parser.parse_args()

#elastic(0) or plastic(1)?
if args.pla:
   pla=1
   print('plasticity was enabled')
else:
   pla=0
   print('plasticity was not enabled')
#cohesive zone width in m: make sure that it converged
#how to calculate the cohesive zone width? 
#run the simulation with fault output for Vr, DS and RF with refining mesh size until the minimum cohesive zone width does not change anymore (very important!)
#formula for the cohesive zone width: (RF-DS)*Vr 
#coh = 760.0
coh = args.coh
#mesh element size
dx = args.dx
#order of accuracy (what you use to compile SeisSol)
order = args.order

print('your minimum cohesive zone is resolved by %f mesh elements' %(coh/dx))
print('your minimum cohesive zone is resolved by %f GPs' %(coh/dx*(order+1)))

#set resolution in context with convergence tests (
# minimum cohesive zone width of 162m (elastic) and 325 (plastic) in convergence tests

if pla==0:
	mincoh = 162.0
else:
	mincoh = 325.0
reso = (dx*mincoh)/coh

#Error slopes for each order, values
#order is always RA, PSR, PSR time, Slip
if pla==0:
	if order==3:
		polycoef_0 = [1.52739978035, 0.836379291573, 1.52201807728, 1.0161358622]
		polycoef_1 = [-4.4893475817, -0.807389129936, -4.34146465619, -2.55381103136]
	elif order==4:
		polycoef_0 = [1.37812937903, 1.11589493374, 1.43176201425, 1.0066205251]
		polycoef_1 = [-4.37561371471, -1.74133399291, -4.41381627321, -2.59964604604]
	elif order==5:
		polycoef_0 = [1.25759353097, 1.32131510275, 1.26380215821, 0.984644208844]
		polycoef_1 = [-4.24036674538, -2.45342685486, -4.11609090554, -2.66217236352]
	elif order==6:
		polycoef_0 = [1.16913715515, 1.38930610536, 1.27245261357, 0.97950635571]
		polycoef_1 = [-4.1179187664, -2.77826477656, -4.17914044254, -2.66901809392]
	else:
		print ('order not defined')

elif pla==1:
	if order==3:
		polycoef_0 = [1.56549445278, 0.726433487514, 1.69615211119, 1.30999916248]
		polycoef_1 = [-4.38905945469, -1.12663481476, -4.68828428197, -3.15076508398]
	elif order==4:
		polycoef_0 = [1.16043927634, 0.848993370158, 1.26072361353, 1.07449236987]
		polycoef_1 = [-3.66553273612, -1.51365000226, -3.89398684582, -2.69926711303]
	elif order==5:
		polycoef_0 = [1.31388913578, 0.916320202014, 1.03772819833, 1.02596911125]
		polycoef_1 = [-4.26237995927, -1.78297284948, -3.5086539198, -2.71310852081]
	elif order==6:
		polycoef_0 = [1.28040171235, 0.932860048315, 1.04397736811, 1.00875346409]
		polycoef_1 = [-4.27329544667, -1.90736538424, -3.5861340745, -2.69759482961]
	else:
		print ('order not defined')

#calculate
error_RA = 10**( polycoef_0[0]*np.log10(reso)+polycoef_1[0] )
error_PSR = 10**( polycoef_0[1]*np.log10(reso)+polycoef_1[1] )
error_PSR_time = 10**( polycoef_0[2]*np.log10(reso)+polycoef_1[2] )
error_Slip = 10**( polycoef_0[3]*np.log10(reso)+polycoef_1[3] )

#output
print('Expected average error (in %) is:')
print('Rupture arrival ', error_RA)
print('Peak slip rate ', error_PSR)
print('Peak slip rate time ', error_PSR_time)
print('Final slip ', error_Slip)
print('RA below 0.2%, PSR below 7% and Final slip below 1% are sufficiently small errors after Day et al. 2005')
