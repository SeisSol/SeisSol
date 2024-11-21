import numpy as np
import os
import argparse

"""
Convert param file to srf file
"""


# srf file requires vs and rho, we take the value from PREM file
PREM_STR = """
  depth radius   vp    vs    rho     Qu   Qk       P
   0.0  6371.0   1.45  0.00  1.02    0.0  57823.0    0.0
   3.0  6368.0   5.80  3.20  2.60  600.0  57823.0    0.0
  15.0  6356.0   6.80  3.90  2.90  600.0  57823.0    0.3
  24.4  6346.6   6.80  3.90  2.90  600.0  57823.0    0.6
  71.0  6300.0   8.08  4.47  3.38  600.0  57823.0    2.2
  80.0  6291.0   8.08  4.47  3.37  600.0  57823.0    2.5
 171.0  6200.0   8.02  4.44  3.36   80.0  57823.0    5.5
 """
if os.path.exists("PREM.txt") is False:
    f = open("PREM.txt", "w")
    f.write(PREM_STR)
    f.close()

PREM = np.loadtxt('PREM.txt', skiprows=1, dtype=float)


def asymmetric_cosine(trise, tfall=None, npts=10000, dt=0.1):
    """
    Initialize a source time function with asymmetric cosine, normalized to 1

    :param trise: rise time
    :type trise: float
    :param tfall: fall time, defaults to trise
    :type trise: float, optional
    :param npts: number of samples
    :type npts: int, optional
    :param dt: sample interval
    :type dt: float, optional
    """
    # initialize
    if not tfall:
        tfall = trise
    t = np.linspace(0, npts * dt, npts, endpoint=False)
    asc = np.zeros(npts)

    # build slices
    slrise = t <= trise
    slfall = np.logical_and(t > trise, t <= trise + tfall)

    # compute stf
    asc[slrise] = 1.0 - np.cos(np.pi * t[slrise] / trise)
    asc[slfall] = 1.0 - np.cos(np.pi * (t[slfall] - trise + tfall) / tfall)

    # normalize
    asc /= trise + tfall

    return asc


def from_usgs_param_file(file, npts, dt, trise_min=0):
    import re
    import pandas as pd
    """
    Reading USGS param file
    Args:
        fh: THE param file
        npts: source time function npts
        dt: dt of source time function
        trise_min: min trise time
        
    Returns:
        nseg (int): no. of segments
        seg_info (list): list of (nx, ny, dx, dy) for each fault segment
        sources (list): list of dataframe of (Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo) for each segments
        stf (list): list of ndarrays of source time function of each subfault in the segments (normalized stf)
    """
    header = "lat lon depth[m] slip[m] rake strike dip t_rup t_ris t_fal mo[Nm]".split()
    fh = open(file, 'r')
    
    lines = fh.readlines()
    if not lines[0].startswith("#Total number of fault_segments"):
        raise ValueError("Not a valid USGS param file.")
    
    
    nseg = int(lines[0].split()[-1])     # number of fault segments
    
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
        seg_info.append([nx, ny, dx, dy])
        
        stf_segments = np.zeros((int(nx*ny), npts))

        line1 = int(istart + 9)
        line2 = int(line1 + nx*ny )
        istart = int(line1 + nx*ny)
        
        # Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo
        data = np.loadtxt(lines[line1:line2])
    
        data[:,2] *= 1e3  # dep: km > m
        data[:,3] *= 1e-2  # slip: cm > m
        data[:,-1] *= 1e-7  # M0: dyn / cm > N * m
    
    
        if np.sum(data[:, 7]< 0) > 0:  # tinit > 0
            raise ValueError(
                "File contains a negative rupture time (unphysical)"
            )
            
        endtime = data[:, 8] + data[:, 8] # t_ris + t_fal
        if np.sum(endtime > (npts - 1) * dt) > 0:
                raise ValueError(
                    "Rise + fall time are longer than the "
                    "total length of calculated slip. "
                    "Please use more samples."
                )

        for j in range(int(nx*ny)):
            trise, tfall = data[j,8], data[j,9]
            _stf = asymmetric_cosine(trise, tfall, npts, dt)
            stf_segments[j,:] = _stf
            
        df = pd.DataFrame(data=data, columns=header)
        sources.append(df)
        stf_all.append(stf_segments)
        
    return nseg, seg_info, sources, stf_all




def find_val(array, val):
    min_val = array[array-val>0].min()
    return min_val


def write_to_srf_file(nseg, seg_info, soruces, stf_all, dt, npts, split_seg=True, srf_file='finite_fault', model_description=""):
    """Write each fault segment to srf file"""
    for iseg in range(nseg):
        df = soruces[iseg]
        no_subfaults = len(df)

        srf = open(f"{srf_file}_seg{iseg+1}.srf",'w')
        srf.write('2.0\n')
        srf.write(f"# {model_description}\n")
        srf.write(f"# nx={seg_info[iseg][0]:.0f} ny={seg_info[iseg][1]:.0f} dx={seg_info[iseg][2]:.2f} dy={seg_info[iseg][3]:.2f}\n")
        srf.write('POINTS %i\n' % no_subfaults)
        
        
        for i in range(no_subfaults):
            LAT, LON, DEP = df.lat.iloc[i], df.lon.iloc[i], df["depth[m]"][i]/1e3
            rake = df.rake.iloc[i]
            rake_rad = np.deg2rad(rake)
            ss_slip, ds_slip = df["slip[m]"].iloc[i] * np.cos(rake_rad)*100, df["slip[m]"].iloc[i] * np.sin(rake_rad)*100 # in cm
            DIP, STRIKE = df.dip.iloc[i], df.strike.iloc[i]
            AREA = seg_info[iseg][2] * seg_info[iseg][2] * 1000*1000*100*100            # in cm^2
            TINIT = df.t_rup.iloc[i]
            DT = dt
            
            ### here we just find the VS and Density from prem model
            ref_depth = find_val(PREM[:,0],DEP)
            Vs = PREM[PREM[:,0]==ref_depth,3][0]  #m/s
            DEN = PREM[PREM[:,0]==ref_depth,4][0]    
            
            subfault_info1 = [LON,LAT,DEP,STRIKE,DIP,AREA,TINIT,DT,Vs*100,DEN]
            subfault_info2 = [rake, ss_slip, npts, ds_slip, npts, 0, 0]
            srf.write('%.3f %.3f %.3f %.1f %.1f %.3e %.3f %.1f %.3e %.3f\n' % tuple(subfault_info1) )
            srf.write('%i %.2f %i  %.2f %i %.2f %i\n' % tuple(subfault_info2))
            _srf = " ".join(f"{x:.2f}" for x in stf_all[iseg][0,:]) ### normalized srf
            srf.write(_srf+"\n")
            srf.write(_srf+"\n")

        srf.close()
        print('Done')
        
        
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Convert USGS file to SRF (2.0)')
    parser.add_argument('filename',type=str, help='USGS param file')  
    parser.add_argument('srf_filename',type=str , help='output file name')
    parser.add_argument('--dt', default=0.5, help='dt of the source time function')  
    parser.add_argument('--npts', default=100, help='npts of the source time function')  
    parser.add_argument('--model_description',default="", help='description of the model stored in SRF file')
    args = parser.parse_args()


    nseg, seg_info, sources, STF = from_usgs_param_file(args.filename, args.npts, args.dt)
    print(f'No. of fault segments in param file: {nseg}')
    write_to_srf_file(nseg, seg_info, sources, STF, args.dt, args.npts, args.model_description, srf_file=args.srf_filename)
