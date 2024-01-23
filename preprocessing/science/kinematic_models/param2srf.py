#!/usr/bin/env python3
import numpy as np
import os
import argparse

"""
Convert param file to srf file
"""


def asymmetric_cosine(trise, tfall=None, dt=0.1):
    """
    Initialize a source time function with asymmetric cosine, normalized to 1

    :param trise: rise time
    :type trise: float
    :param tfall: fall time, defaults to trise
    :type trise: float, optional
    :param dt: sample interval
    :type dt: float, optional
    """
    # initialize
    if not tfall:
        tfall = trise
    t = np.linspace(0, trise+tfall, int((trise+tfall)/ dt)+1, endpoint=True)
    asc = np.zeros_like(t)

    # build slices
    slrise = t <= trise
    slfall = np.logical_and(t > trise, t <= trise + tfall)

    # compute stf
    asc[slrise] = 1.0 - np.cos(np.pi * t[slrise] / trise)
    asc[slfall] = 1.0 - np.cos(np.pi * (t[slfall] - trise + tfall) / tfall)

    # normalize
    integral_aSTF = np.trapz(np.abs(asc[:]), dx=dt)
    asc /= integral_aSTF

    return asc


def from_usgs_param_file(file, dt, trise_min=0):
    import re
    import pandas as pd
    from io import StringIO
    """
    Reading USGS param file
    Args:
        fh: THE param file
        dt: dt of source time function
        trise_min: min trise time
        
    Returns:
        nseg (int): no. of segments
        seg_info (list): list of (nx, ny, dx, dy) for each fault segment
        sources (list): list of dataframe of (Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo) for each segments
        stf (list): list of ndarrays of source time function of each subfault in the segments (normalized stf)
    """
    header = "lat lon depth slip rake strike dip t_rup t_ris t_fal mo"
    with open(file, 'r') as fid:
        lines = fid.readlines()

    if not lines[0].startswith("#Total number of fault_segments"):
        raise ValueError("Not a valid USGS param file.")
    
    
    nseg = int(lines[0].split()[-1])     # number of fault segments
    print(f'No. of fault segments in param file: {nseg}')
    
    seg_info = []
    sources = []
    stf_all = []


    fault_seg_line = [l for l in lines if "#Fault_segment " in l]
    assert len(fault_seg_line)==nseg, f"No. of segments are wrong. {len(fault_seg_line)} {nseg}"
    
    istart = 1
    for i_seg in range(nseg):
        numbers = re.findall(r'[-+]?\d*\.\d+|\d+', fault_seg_line[i_seg])

        # Convert extracted strings to float
        numbers = [float(num) for num in numbers]
        _, nx, dx, ny, dy = numbers
        nx, ny = int(nx), int(ny)

        seg_info.append([nx, ny, dx, dy])       
        stf_segments = []

        line1 = istart + 9
        line2 = line1 + nx*ny
        istart = line1 + nx*ny

        text_file = StringIO("\n".join([header, *lines[line1:line2]]))   
        df = pd.read_csv(text_file, sep='\s+')
        
        assert (df['t_rup'] >= 0).all(), "AssertionError: Not all rupture time are greater than or equal to 0."

        for j in range(int(nx*ny)):
            _stf = asymmetric_cosine(df['t_ris'][j], df['t_fal'][j], dt)
            stf_segments.append(_stf)
            
        sources.append(df)
        stf_all.append(stf_segments)
        
    return seg_info, sources, stf_all


def write_to_srf_file(seg_info, sources, stf_all, dt, split_seg=True, srf_file='finite_fault', model_description=""):
    """Write each fault segment to srf file"""
    nseg = len(seg_info)
    for iseg in range(nseg):
        df = sources[iseg]
        no_subfaults = len(df)

        top_center_lon, top_center_lat = df.lon.iloc[0], df.lat.iloc[0]
        nx, ny, dx, dy = seg_info[iseg]
        lx = nx * dx
        ly = ny * dy
        strike = df.strike.iloc[0] 
        dip = df.dip.iloc[0]
        index_hypo = df['t_rup'].idxmin()
        j_hypo, i_hypo = divmod(index_hypo, nx)
        along_strike_location_hypo = i_hypo * dx
        along_dip_location_hypo = j_hypo * dy

        with open(f"{srf_file}_seg{iseg+1}.srf",'w') as srf:
            srf.write('2.0\n')
            srf.write("PLANE 1\n")
            srf.write(f"{top_center_lon} {top_center_lat} {nx} {ny} {lx} {ly} ")
            srf.write(f"{strike} {dip} {along_strike_location_hypo} {along_dip_location_hypo}\n")
            srf.write(f'POINTS {no_subfaults}\n')
            
            for i in range(no_subfaults):
                LAT, LON, DEP = df.lat.iloc[i], df.lon.iloc[i], df["depth"][i]
                rake = df.rake.iloc[i]
                rake_rad = np.deg2rad(rake)
                DIP, STRIKE = df.dip.iloc[i], df.strike.iloc[i]
                AREA = seg_info[iseg][2] * seg_info[iseg][2] * 1000*1000*100*100            # in cm^2
                TINIT = df.t_rup.iloc[i]
                DT = dt
                
                Vs = 0.0
                DEN = 0.0
                subfault_info1 = [LON,LAT,DEP,STRIKE,DIP,AREA,TINIT,DT,Vs,DEN]
                subfault_info2 = [rake, df.slip.iloc[i], len(stf_all[iseg][i]), 0, 0, 0, 0]
                srf.write('%.3f %.3f %.3f %.1f %.1f %.3e %.3f %.1f %.3e %.3f\n' % tuple(subfault_info1) )
                srf.write('%i %.2f %i  %.2f %i %.2f %i\n' % tuple(subfault_info2))
                stf = stf_all[iseg][i]
                #format string into groups of 10 values for readability
                _srf = "\n".join(" ".join(f"{x:.5f}" for x in stf[i:i+10]) for i in range(0, len(stf), 10))
                srf.write(_srf+"\n")

        print('Done')
        
        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Convert USGS file to SRF (2.0)')
    parser.add_argument('filename',type=str, help='USGS param file')  
    parser.add_argument('srf_filename',type=str , help='output file name')
    parser.add_argument('--dt', default=0.5, help='dt of the source time function')  
    parser.add_argument('--model_description',default="", help='description of the model stored in SRF file')
    args = parser.parse_args()


    seg_info, sources, STF = from_usgs_param_file(args.filename, args.dt)
    write_to_srf_file(seg_info, sources, STF, args.dt, args.model_description, srf_file=args.srf_filename)
