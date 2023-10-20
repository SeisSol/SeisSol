#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 21:18:55 2023

@author: fkutschera
"""
import matplotlib
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import os
import glob

description = '''Plot SeisSol energy output(s). 
You can run the script with the default (Moment rate [MR]) or select multiple energy outputs available from a SeisSol simulation.
You can also compare the energy outputs of different simulations.)'''
parser = argparse.ArgumentParser(description=description)
parser.add_argument("Input", nargs="+", metavar=("file"),
                    help="energy.csv file(s) from your SeisSol simulation. Multiple input files possible.")
parser.add_argument("-d", "--Data", type=str, required=False, nargs="+", metavar=("variable"), default=["MR"], 
                    help="Data to plot. Default: Moment rate (MR). Other available options are Gravitational energy (GE), Acoustic energy (AE), Acoustic kinetic energy (AKE), Elastic kinetic energy (EKE), Elastic energy (EE), Total frictional work (TFW), Static frictional work (SFW), Seismic moment (SM), Plastic moment (PM)")
parser.add_argument("--xlim", nargs="+", type=float, help="xmin xmax")
parser.add_argument("--ylim", nargs="+", type=float, help="ymin ymax")
parser.add_argument("--xlabel", default="Time [s]", type=str, help="Choose xlabel")
parser.add_argument("--ylabel", default="", type=str, help="Choose ylabel")
parser.add_argument("--dpi", default=300, type=int, help="Dots per Inch")
parser.add_argument("-o", "--Output", help='Name of the output file', default='energy_output.png')
args = parser.parse_args()


### search through wildcards
files = [] 
for f in args.Input:
    files.extend(glob.glob(f))
    
for i in range(0,len(files)):
    if os.path.isfile(os.path.expanduser(args.Input[i])):
        print("File selected: {}".format(os.path.expanduser(args.Input[i])))
    else:
        sys.exit("File '{}' does not exist.".format(os.path.expanduser(args.Input[i])))
        
    df = pd.read_csv(args.Input[i])
    M0 = df['seismic_moment'][len(df)-1]
    Mw = 2.0/3.0 * np.log10(M0) - 6.07;
    # see https://github.com/SeisSol/SeisSol/blob/829f5e258c12046ae3626837bab6d19d13d6c006/src/ResultWriter/EnergyOutput.cpp#L436C45-L436C45

    print("M0 = {:.2e} and Mw = {:.2f}".format(M0, Mw))
        
    for j in range(0,len(args.Data)):
        if (args.Data[j].upper() == "MR"):
            energy_output = "moment_rate"
            df[energy_output] = np.gradient(df["seismic_moment"], df["time"].iloc[1])
        elif (args.Data[j].upper() == "GE"):
            energy_output = "gravitational_energy"
        elif (args.Data[j].upper() == "AE"):
            energy_output = "acoustic_energy"
        elif (args.Data[j].upper() == "AKE"):
            energy_output = "acoustic_kinetic_energy"
        elif (args.Data[j].upper() == "EE"):
            energy_output = "elastic_energy"
        elif (args.Data[j].upper() == "EKE"):
            energy_output = "elastic_kinetic_energy"
        elif (args.Data[j].upper() == "TFW"):
            energy_output = "total_frictional_work"
        elif (args.Data[j].upper() == "SFW"):
            energy_output = "static_frictional_work"
        elif (args.Data[j].upper() == "SM"):
            energy_output = "seismic_moment"
        elif (args.Data[j].upper() == "PM"):
            energy_output = "plastic_moment"
        else:
            sys.exit("Are you sure? Please select minimum one of the following: MR, GE, AE, AKE, EKE, EE, TFW, SFW, SM, PM")

        print("Energy output selected is {}: {}".format(args.Data[j].upper(), energy_output))
        
        plt.plot(df["time"], df[energy_output], label="Mw {:.2f}, {}".format(Mw,args.Data[j].upper()))
    
        if (args.ylabel != ""):
            ylabel = args.ylabel
        elif (args.Data[j].upper() == "MR") and (int(len(args.Data)) == 1):
            ylabel = "Moment release [Nm s$^{-1}$]"
        elif (args.Data[j].upper() == "PM") and (int(len(args.Data)) == 1):
            ylabel = "Plastic moment [Nm]"
        elif (args.Data[j].upper() == "SM") and (int(len(args.Data)) == 1):
            ylabel = "Seismic moment [Nm]"
        elif ("energy" in energy_output or "work" in energy_output):
            ylabel = "Energy [J]"
        else:
            ylabel = args.ylabel
        
    plt.xlabel(args.xlabel)
    plt.ylabel(ylabel)  

    plt.legend()
        
    if (args.xlim != None):
        plt.xlim(args.xlim[0], args.xlim[1])
    if (args.ylim != None):
        plt.ylim(args.ylim[0], args.ylim[1])
        
plt.savefig("{}".format(args.Output), dpi=args.dpi)
