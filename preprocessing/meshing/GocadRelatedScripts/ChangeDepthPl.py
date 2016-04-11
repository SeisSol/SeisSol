#Author: Thomas Ulrich, LMU
 
# parsing python arguments
import argparse
import os
parser = argparse.ArgumentParser(description='read a pl file and change its depth to a constant value')
parser.add_argument('input_file', help='.pl file with 2 close curves')
parser.add_argument('output_file', help='gocad .pl output file with desired depth')
parser.add_argument('--depth', nargs=1, metavar=('d'), default = ('1000'), help='desired depth (absolute value)')
args = parser.parse_args()


fid = open(args.input_file)
lines = fid.readlines()
fid.close()

fout = open(args.output_file,'w')

for line in lines:
   if line.startswith('VRTX'):
      listsplit = line.split()[0:-1]
      listsplit.append('-'+args.depth[0])
      mystring = "%s %s %s %s %s\n" % tuple(listsplit)
      fout.write(mystring)
   else:
      fout.write(line)
fout.close()
