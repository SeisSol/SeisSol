import os.path

myfolder = 'OutputParRF10-25-100/'
myprefix = 'data'
PrintInterval=1500
dt=0.000306614
tf=20.

fout = open('newpvd.pvd','w')

fout.write("<?xml version='1.0' ?>\n<VTKFile type='Collection' version='0.1'>\n<Collection>\n")


iPrint=1

while(True):
	print(dt*iPrint)
	iPart=0
	for icore in range(0,256):
		ti=iPrint*dt
		fname='%s-fault-00%03d' %(myprefix, icore)
		if os.path.isdir(myfolder + fname):
			fout.write("<DataSet timestep='%f' group='' part='%d' file='%s-fault-00%03d/timestep-00000%04d.vtu'/>\n" %(ti, iPart, myprefix, icore, iPrint))
			iPart = iPart+1
	if (iPrint== (int) (tf/dt)):
		break
	iPrint = iPrint + PrintInterval
	if (iPrint>tf/dt):
		iPrint=(int) (tf/dt)

fout.write("</Collection>\n</VTKFile>\n")


fout.close()
	

