#convert inp (from SimModeler 2d mesh (ABAQUS 2D)) to ts (Gocad)

# parsing python arguments
import argparse
parser = argparse.ArgumentParser(description='convert inp (from SimModeler5 2d mesh (ABAQUS 2D)) to ts (Gocad)')
parser.add_argument('inp_filename', help='inp filename (SimModeler5 2d mesh (ABAQUS 2D))')
parser.add_argument('ts_filename', nargs='?',help='output filname (if not used = inpbasename.ts)',default='')
parser.add_argument('--isolate', dest='isolate', action='store_true', help='isolate every Gocad surface in a different stl solid')
args = parser.parse_args()

if args.ts_filename == '':
   args.ts_filename = args.inp_filename[0:-4]+'.ts'

fid = open(args.inp_filename)
lines = fid.readlines()
fid.close()

for i,line in enumerate(lines):
   if line.startswith('*Node'):
      inodes=i+1
   if line.startswith('*Element'):
      iel=i+1
      break

fout = open (args.ts_filename, 'w')

mystring = "GOCAD TSURF 1\n\
HEADER {\n\
name:s1\n\
}\n\
TRIANGLES\n"

for i in range(inodes, iel-1):
   val = lines[i].split(',')
   mystring = mystring + 'VRTX %s %s %s %s' %(val[0], val[1],val[2],val[3])

fout.write(mystring)

for i in range(iel, len(lines)):
   line = lines[i]
   if line.startswith('*Element'):
      if args.isolate:
         fout.write('END\n'+mystring)
      continue 
   val = lines[i].split(',')
   fout.write('TRGL %s %s %s' %(val[1],val[2],val[3]))
fout.write('END')
fout.close()
