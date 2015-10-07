import numpy as np
import math
import os

#read Gambit mesh

fid=open('tpv31c2c105-15-500-5000.neu')

for i in range(0,6):
   fid.readline()

l7=fid.readline()
data=[int(el) for el in l7.split()]
NUMNP = data[0]
NELEM = data[1]
NGRPS = data[2]
NBSETS = data[3]
print("%d nodes %d elements %d groups %d sets" %(NUMNP, NELEM, NGRPS, NBSETS))

fid.readline()
fid.readline()

for i in range(0,NUMNP):
   fid.readline()
print("nodes read")

fid.readline()
fid.readline()

elems=np.zeros((NELEM+1,4))
for i in range(1,NELEM+1):
   data=fid.readline()
   elems[i]=[int(el) for el in data.split()[3:]]
print(elems)

fmesh=open("elems.dat",'w')
fmesh.write("%d\n" %NELEM)
np.savetxt(fmesh,elems[1:], fmt='%d %d %d %d')
fmesh.close()

#creating metis graph to access more easily the neigbhors
print("running m2gmetis")
cmd='m2gmetis elems.dat elems.graph -ncommon=3'
os.system( cmd )
print("done")

#load the graph in a variable
fg=open("elems.graph")
lines=fg.readlines()
fg.close()
graph=np.zeros((NELEM+1,4))

for i in range(1,NELEM+1):
   line=lines[i]
   data=[int(el) for el in line.split()]
   for j in range(len(data)):
       graph[i,j]=data[j]
print graph

fid.readline()
fid.readline()

#read groups
for i in range(0,NGRPS):
   nel= int(fid.readline().split()[3])
   print("reading group %d composed of %d elements" %(i,nel))
   fid.readline()
   fid.readline()
   for i in range(0,int(math.ceil(nel/10.))):
      fid.readline()

#boundary 101
print(fid.readline())
fid.readline()

data=fid.readline().split()
boundtype = int(data[0])
nel = int(data[2])
print("processing boundary %d composed of %d elements" % (boundtype, nel))
for i in range(0,nel):
   fid.readline()

#boundary 103
fid.readline()
fid.readline()

data=fid.readline().split()
boundtype = int(data[0])
nel = int(data[2])
print("processing boundary %d composed of %d elements" % (boundtype, nel))
bound103=np.zeros((nel,2))

for i in range(0,nel):
   data=[int(el) for el in fid.readline().split()]
   bound103[i]=[data[0],data[2]]
   print(bound103[i])

   nodes=elems[bound103[i,0]]
   faceid=bound103[i,1]
   neighElems = graph[[bound103[i,0]]]

   if faceid == 1:
      facenodes=[nodes[0], nodes[2], nodes[1]]
   elif faceid == 2:
      facenodes=[nodes[0], nodes[1], nodes[3]]
   elif faceid == 3:
      facenodes=[nodes[0], nodes[3], nodes[2]]
   elif faceid == 4:
      facenodes=[nodes[1], nodes[2], nodes[3]]
   else:
      print("unknown faceid %d" %faceid)
      stop
   for neiEl in neighElems:
      if el==0:
         break
      el=elems[neiEl]
      if (facenodes[0] in el) & (facenodes[1] in el) & (facenodes[2] in el):
         if not (el[0] in facenodes):
            faceidel=4
         elif not (el[1] in facenodes):
            faceidel=3
         elif not (el[2] in facenodes):
            faceidel=2
         elif not (el[3] in facenodes):
            faceidel=1
         else:
            print("unable de determine faceid")
            stop
         print("%d %d " %(i+1,faceidel))
          
print(bound103)





#boundary 105
fid.readline()
fid.readline()

data=fid.readline().split()
boundtype = int(data[0])
nel = int(data[2])
print("processing boundary %d composed of %d elements" % (boundtype, nel))
for i in range(0,nel):
   fid.readline()

