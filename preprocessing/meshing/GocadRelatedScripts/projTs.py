import argparse
parser = argparse.ArgumentParser(description='project ts file')
parser.add_argument('ts_file', help='ts filename')
parser.add_argument('--proj', nargs=1, metavar=('projname'), default = (''), help='string describing its projection (ex: +init=EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered')
args = parser.parse_args()

#set projection
import mpl_toolkits.basemap.pyproj as pyproj
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
if args.proj[0]!='geocent':
   sProj = args.proj[0]
   myproj=pyproj.Proj(sProj)
else:
   myproj = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')

#read Ts file
fid=open(args.ts_file)
lines = fid.readlines()
fid.close()

foutname = args.ts_file[0:-3]+'_proj.ts'
fout=open(foutname,'w')

for line in lines:
   if line.startswith('VRTX'):
      val=[float(val) for val in line.split()[1:5]]
      xyz = pyproj.transform(lla, myproj, val[1],val[2], val[3]*1e3, radians=False)
      fout.write('VRTX '+str(int(val[0]))+ ' %.10e %.10e %.10e\n' %tuple(xyz))
   else:
      fout.write(line)
fout.close()

