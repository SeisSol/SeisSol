# Script to determine the rise time (time the absolute slip rate is above a certain threshold) from seissol fault.xdmf output
# Output of this script is a new xdmf file

import numpy as np
import seissolxdmf
import argparse
import os

parser = argparse.ArgumentParser(description='calculate rise times from fault xdmf file')
parser.add_argument('filename', help='path+filename-fault.xdmf')
parser.add_argument('--threshold', nargs=1, metavar=('threshold'), default=([0.01]), help='specify slip rate threshold' ,type=float)
parser.add_argument('--dt', nargs=1, metavar=('dt'), default=([0.1]), help='insert time step of xdmf file' ,type=float)
parser.add_argument('--events', nargs=1, metavar=('events'), choices=[1,2], default=([1]), help='number of events [1,2] in xdmf file' ,type=int)
args = parser.parse_args()

def counting_loop(row, array, endcol, startcol=0, threshold=0.01):
    count = 0
    for j in range(startcol, endcol):
        if array[row,j] < threshold and count > 0: break
        if array[row,j] <  threshold and count == 0: continue
        if array[row,j] >= threshold: count += 1
    return count

def get_rise_time(ASR_array, threshold=0.01, dt=0.1, events=1):
    if events == 1:
        RT = np.zeros(ASR_array.shape[1])
        ASR_array=np.transpose(ASR_array) # to speed up
        for i in range(ASR_array.shape[0]):
            count = counting_loop(i, ASR_array, ASR_array.shape[1], threshold=threshold) 
            RT[i] = (count * dt)
            
    if events == 2:
        RT = np.zeros_like(ASR_array[0:2,:])
        mid = int(ASR_array.shape[0]/2)
        ASR_array=np.transpose(ASR_array)
        for i in range(ASR_array.shape[0]):
            count = counting_loop(i, ASR_array, mid, threshold=threshold)      
            RT[0,i] = (count * dt) 
            count = counting_loop(i, ASR_array, ASR_array.shape[1], startcol=mid, threshold=threshold)    
            RT[1,i] = (count * dt)
            
    return RT

def get_absolute_slip_rate(sx_object):
    SRs = sx_object.ReadData('SRs')
    SRd = sx_object.ReadData('SRd')
    ASR = np.sqrt(SRs**2 + SRd**2)
    return ASR

def CreateXdmf(fn, nNodes, ntriangles, aDataName):     
    # original by Thomas Ulrich:
    # https://github.com/SeisSol/Meshing/blob/master/vizualizeBoundaryConditions/vizualizeBoundaryConditions.py
    xdmf="""<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
 <Domain>
  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">
   <Grid Name="step_000000000000" GridType="Uniform"><!-- mesh id: 0, mesh step: 0 -->
    <Topology TopologyType="Triangle" NumberOfElements="%d">
     <DataItem NumberType="Int" Precision="8" Format="HDF" Dimensions="%d 3">%s.h5:/mesh0/connect</DataItem>
    </Topology>
    <Geometry name="geo" GeometryType="XYZ" NumberOfElements="%d">
     <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="%d 3">%s.h5:/mesh0/geometry</DataItem>
    </Geometry>
    <Time Value="0"/>""" %(ntriangles, ntriangles, fn, nNodes, nNodes, fn)

    for dataName in aDataName:
        xdmf=xdmf + """
    <Attribute Name="%s" Center="Cell" Type="Scalar">
      <DataItem NumberType="Float" Precision="8" Format="HDF" Dimensions="1 %d">%s.h5:/mesh0/%s
     </DataItem>
    </Attribute>""" % (dataName,ntriangles,fn,dataName)

    xdmf=xdmf + """
   </Grid>
  </Grid>
 </Domain>
</Xdmf>
"""
    fid=open(fn+'.xdmf','w')
    fid.write(xdmf)
    fid.close()
    
def write_xdmfh5(fname, aDataName, xyz, connect, RT):
    # original by Thomas Ulrich:
    # https://github.com/SeisSol/Meshing/blob/master/vizualizeBoundaryConditions/vizualizeBoundaryConditions.py
    nNodes=xyz.shape[0]
    ntriangles = connect.shape[0]
    CreateXdmf(fname, nNodes, ntriangles, aDataName)
    #Write h5 file
    import h5py
    h5f = h5py.File(fname+'.h5','w')
    h5f.create_dataset('mesh0/connect', (ntriangles,3), dtype='i8')
    h5f['mesh0/connect'][:,:] = connect[:,:]
    h5f.create_dataset('mesh0/geometry', xyz.shape, dtype='d')
    h5f['mesh0/geometry'][:,:] = xyz[:,:]
    for i in range(0, len(aDataName)):
        hdname = "mesh0/"+aDataName[i]
        h5f.create_dataset(hdname, (1, ntriangles), dtype='d')
        h5f[hdname][0,:] = eval(aDataName[i])[:]
    h5f.close()
    print ("done writing %s.xdmf" %fname)

print("Loading file and calculating absolute slip rates...") 
sx = seissolxdmf.seissolxdmf(args.filename)
ASR = get_absolute_slip_rate(sx)
          
print("Determining rise time...")
RT = get_rise_time(ASR, threshold=args.threshold[0], dt=args.dt[0], events=args.events[0])

print("Writing output file...")
connect = sx.ReadConnect()
xyz = sx.ReadGeometry()
if args.events[0] == 1:
    aDataName = ['RT']
else:
    aDataName = ['RT1', 'RT2']
    RT1 = RT[0,:]
    RT2 = RT[1,:]
prefix, ext = os.path.splitext(args.filename)
fn = prefix+'_RT'
write_xdmfh5(fn, aDataName, xyz, connect, RT)
