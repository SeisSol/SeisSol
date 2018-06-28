import numpy as np
import argparse

parser = argparse.ArgumentParser(description='project ts file')
parser.add_argument('ts_file', help='ts filename')
parser.add_argument('--proj', nargs=1, metavar=('projname'), default = (''), help='string describing its projection (ex: +init=EPSG:32646 (UTM46N), or geocent (cartesian global)) if a projection is considered')
parser.add_argument('--debug', dest='debug', action='store_true', help='print additionnal debug info')
parser.add_argument('--merged', dest='merged', action='store_true', help='if true, each surface is not isolated in a stl solid')
parser.add_argument('--bstl', dest='bstl', action='store_true', help='if true, write binary stl instead of stl (note that surface will be merged once imported by SimModeler)')
args = parser.parse_args()


#set projection
if args.proj!='':
   import mpl_toolkits.basemap.pyproj as pyproj
   lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
   if args.proj[0]!='geocent':
      sProj = args.proj[0]
      myproj=pyproj.Proj(sProj)
   else:
      myproj = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')

#set logging
import logging, sys
if args.debug:
   logging.basicConfig(format='%(message)s',stream=sys.stderr, level=logging.DEBUG)
else:
   logging.basicConfig(format='%(message)s',stream=sys.stderr, level=logging.INFO)

surface_name='undefined'
sid=0

if not args.bstl:
   foutname = args.ts_file[0:-3]+'.stl'
   fout = open(foutname,'w')
   if args.merged:
      fout.write("solid 0\n")
else:
   import struct
   foutname = args.ts_file[0:-3]+'.bstl'
   fout = open(foutname,'wb')
   fout.seek(80)
   #number of triangle. will be filled afterwards
   fout.write(struct.pack('<L', 0))


nfacets=0
with open(args.ts_file) as fid:
   while True:
      xyzl=[]
      vid={}
      trl=[]
      prev_vert=-1
      ivertex=0
      line=fid.readline()
      if not line:
         break
      title = line.split()
      #Skip everything not Surface
      if title[1].lower()!='tsurf':
         logging.info("skipping %s" %(title[1]))
         while True:
            line=fid.readline()
            if not line:
               break
            if line.startswith('END_'):
               continue
            elif line.startswith('END'):
               break
      else:
         while True:
            line=fid.readline()
            if not line:
               break
            if line.startswith('VRTX'):
               val=[float(val) for val in line.split()[1:5]]
               newVid = int(val[0])
               vid[newVid] =  ivertex
               ivertex=ivertex+1
               xyzl.append(val[1:4])
            elif line.startswith('TRGL'):
               val=[int(val) for val in line.split()[1:4]]
               trl.append(val)
            elif line.startswith('END_'):
               continue
            elif line.startswith('name'):
               surface_name = (line.split(':')[1]).strip()
               logging.info("now processing %s" %surface_name)
               if not args.merged:
                  sid=sid+1
            elif line.startswith('ATOM'):
               val=[int(val) for val in line.split()[1:3]]
               vid0 = vid[val[1]]
               xyzl.append(xyzl[vid0])
               vid[val[0]] =  ivertex
               ivertex=ivertex+1
            elif line.startswith('END'):
               break
         nodes = np.asarray(xyzl)
         triangles = np.asarray(trl)
         ntriangles = np.shape(triangles)[0]
         logging.debug("reindexing triangles...")
         for itr in range(ntriangles):
            for iv in range(3):
               triangles[itr,iv] = vid[triangles[itr,iv]]

         nnodes=np.shape(nodes)[0]
         nfacets = nfacets + ntriangles
         logging.debug("done reading %s, found %d nodes and %d triangles" %(surface_name, nnodes,ntriangles))
         if args.proj!='':
            logging.info("projecting the nodes coordinates")
            xyzb = pyproj.transform(lla, myproj, nodes[:,0],nodes[:,1], 1e3*nodes[:,2], radians=False)
            nodes[:,0]= xyzb[0]
            nodes[:,1]= xyzb[1]
            nodes[:,2]= xyzb[2]


         #compute efficiently the normals
         logging.debug("computing the normals")
         normal = np.cross(nodes[triangles[:,1],:]-nodes[triangles[:,0],:],nodes[triangles[:,2],:]-nodes[triangles[:,0],:])

         #small change to make it compatible with older version of Python
         #norm=np.linalg.norm(normal, axis=1)
         norm=np.apply_along_axis(np.linalg.norm, 1, normal)
         
         logging.debug("removing flat triangles")
         #This is not necessary but then it generate Warnings in this script for 0 division
         #and in SimModeler import than destabilize new users
         ids_flat_triangles = np.where(norm>0)
         norm = norm[ids_flat_triangles]
         normal = normal[ids_flat_triangles]
         triangles = triangles[ids_flat_triangles]
         ntriangles = np.shape(norm)[0]

         normal = normal/norm.reshape((ntriangles,1))

         logging.debug("Writing solid %s" %surface_name)
         if not args.bstl:
            if not args.merged:
               fout.write("solid %s\n" %surface_name)
            for k in range(ntriangles):
                 fout.write('facet normal %e %e %e\n' %tuple(normal[k,:]))
                 fout.write('outer loop\n')
                 for i in range(0,3):
                    fout.write('vertex %.10e %.10e %.10e\n' % tuple(nodes[triangles[k,i],:]))
                 fout.write("endloop\n")
                 fout.write("endfacet\n")
            if not args.merged:
               fout.write("endsolid %s\n" %surface_name)
         else:
            for k in range(ntriangles):
               fout.write(struct.pack('<3f', *normal[k,:]))
               for i in range(0,3):
                  fout.write(struct.pack('<3f', *nodes[triangles[k,i],:]))
               fout.write(struct.pack('<H', sid))
         logging.debug("done writing the file")
if not args.bstl:
   if args.merged:
      fout.write("endsolid 0\n")
else:
   fout.seek(80)
   fout.write(struct.pack('<L', nfacets))
fout.close()
