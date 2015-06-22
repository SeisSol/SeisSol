#!/usr/bin/python
##
# @file
# This file is part of SeisSol.
#
# @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
#
# @section LICENSE
# Copyright (c) 2013, SeisSol Group
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import collections
import math
import re

FLOAT_REGEX_STRING = '[+\-0-9\.eE]'
ENDSECTION = 'ENDOFSECTION'
ENDSECTION_REGEX = re.compile('^ENDOFSECTION$')

class GambitFile(collections.Iterable):
    """Provides an file like iterator but allows seek/tell and removes comments/empty lines"""
    
    def __init__(self, file):
        self.__file = file
        
    def __del__(self):
        self.__file.close()
        
    def __iter__(self):
        line = self.__file.readline()
        while line:
            line = GambitFile.__prepareLine(line)
            if line == False:
                continue
            
            yield line
            
            line = self.__file.readline()
    
    def seek(self, pos):
        self.__file.seek(pos)
        
    def tell(self):
        return self.__file.tell()
    
    def skipSection(self):
        for line in self:
            # End of section?
            if line == ENDSECTION:
                break
    
    @staticmethod
    def __prepareLine(line):
        # Comment?
        if line[0] == '/':
            return False
            
        return line.strip()

class GambitReader(object):
    """Reads a gambit neutral file"""
    
    class Cells(collections.Iterable):
        def __init__(self, len, file):
            self.__len = len
            self.__file = file
            self.__seek = file.tell()
            
            # Iterator finished?
            self.__finished = True
            
            # Skip cells now, so the reader can continue
            file.skipSection()
                
        def __len__(self):
            return self.__len
        
        def __iter__(self):
            file = self.__file
            
            # Reset iterator
            if self.__finished:
                self.__finished = False
                file.seek(self.__seek)
                
            # This regex only supports tetrahedra
            tetCell = re.compile('^('+FLOAT_REGEX_STRING+'+) +6 +4 +(\d+) +(\d+) +(\d+) +(\d+)$')
                
            nextCellNumber = 1
            for line in file:
                # Check for cell
                match = re.match(tetCell, line)
                if match:
                    # Check for correct cell number
                    if float(match.group(1)) != nextCellNumber:
                        raise IOError('Cells are not sorted, wrong cell number after '+str(nextCellNumber))
                    nextCellNumber += 1
                
                    yield [int(match.group(x))-1 for x in [2, 3, 4, 5]]
                    continue
            
                # End of section?
                if line == ENDSECTION:
                    self.__finished = True
                    return
            
                raise IOError('Error in cell section in line: '+line)
            
    
    def __init__(self, file, sparse = False):
        """sparse currently does not support coords/groups/boundaries"""
        
        # Open file if a string is specified
        if isinstance(file, basestring):
            file = open(file, 'rU')
            
        # We need to close it later
        self.__file = GambitFile(file)
            
        # Read header
        self.__readHeader(self.__file)
        
        # Create array for groups and boundaries
        self.__groups = []
        self.__boundaries = []
        
        # Read section by section
        for line in self.__file:
            if line == 'NODAL COORDINATES '+self.version:
                self.__readCoords(self.__file, not sparse)
            elif line == 'ELEMENTS/CELLS '+self.version:
                self.__readCells(self.__file, sparse)
            elif line == 'ELEMENT GROUP '+self.version:
                self.__readGroup(self.__file, not sparse)
            elif line == 'BOUNDARY CONDITIONS '+self.version:
                self.__readBoundary(self.__file, not sparse)
            else:
                raise IOError('Could not parse line: '+line)
        
        # Check header sizes
        if self.problemSize[0] != len(self.__coords):
            raise IOError('Number of coordinates not correct: number in header = '
                +str(self.problemSize[0])+' != actual number = '+str(len(self.__coords)))

        if self.problemSize[1] != len(self.__cells):
            raise IOError('Number of cells not correct: number in header = '
                +str(self.problemSize[1])+' != actual number = '+str(len(self.__cells)))

        if self.problemSize[2] != len(self.__groups):
            raise IOError('Number of groups not correct: number in header = '
                +str(self.problemSize[2])+' != actual number = '+str(len(self.__groups)))
            
        if self.problemSize[3] != len(self.__boundaries):
            raise IOError('Number of boundaries not correct: number in header = '
                +str(self.problemSize[3])+' != actual number = '+str(len(self.__boundaries)))
            
    def __del__(self):
        del self.__file
            
    def __readHeader(self, file):
        header = collections.deque([
            (re.compile('^CONTROL INFO (\d+\.\d+\.\d+)$'), 'version', lambda m: m.group(1)),
            (re.compile('^\*\* GAMBIT NEUTRAL FILE$'), None),
            (re.compile('^(.*)$'), 'name', lambda m: m.group(1)),
            (re.compile('^PROGRAM:  .{20} +(VERSION: +)?\d+\.\d+\.\d+$'), None),
            (re.compile('^.*$'), 'date', lambda m: m.group(0)),
            (re.compile('^NUMNP {5}NELEM {5}NGRPS {4}NBSETS {5}NDFCD {5}NDFVL$'), None),
            (re.compile('^(\d+) +(\d+) +(\d+) +(\d+) +(\d+) +(\d+)$'),
                'problemSize', lambda m: [int(m.group(x+1)) for x in range(6)]),
            (ENDSECTION_REGEX, None)
        ])
                
        for line in file:
            # Check if we found the correct header
            match = re.match(header[0][0], line)
            if not match:
                raise IOError('Could not parse header line: '+line)
                
            # Save matched expression?
            if header[0][1]:
                self.__setattr__(header[0][1], header[0][2](match))
                
            # Header found -> goto next header
            header.popleft()
                    
            # All headers -> stop
            if not header:
                break
            
    def __readCoords(self, file, save):
        if save:
            coord = re.compile('^('+FLOAT_REGEX_STRING+'+) +('+FLOAT_REGEX_STRING+'+) +('
                +FLOAT_REGEX_STRING+'+) +('+FLOAT_REGEX_STRING+'+)$')
            
            self.__coords = []
            
            nextCoordNumber = 1
            for line in file:
                # Check for a coordinate
                match = re.match(coord, line)
                if match:
                    # Check for correct coordinate number
                    if float(match.group(1)) != nextCoordNumber:
                        raise IOError('Coordinates are not sorted, wrong coordinate number after '+str(nextCoordNumber))
                    nextCoordNumber += 1
                     
                    self.__coords.append([float(match.group(x)) for x in [2, 3, 4]])
                    continue
                  
                # End of section?
                if line == ENDSECTION:
                    break
                
                raise IOError('Error in coordinate section in line: '+line)
        else:
            self.__coords = xrange(self.problemSize[0])
            file.skipSection()
        
    def __readCells(self, file, sparse):
        if sparse:
            self.__cells = GambitReader.Cells(self.problemSize[1], file)
            return
        
        self.__cells = [c for c in GambitReader.Cells(self.problemSize[1], file)]
        
    def __readGroup(self, file, save):
        if save:
            header = collections.deque([
                (re.compile('^GROUP: +\d+ ELEMENTS: +(\d+) ?MATERIAL: +(\d+) NFLAGS: +\d+$'),
                    'header', lambda m: [int(m.group(i)) for i in [1, 2]]),
                (re.compile('^.*$'), 'name', lambda m: m.group(0)),
                (re.compile('^.*$'), None) # Ignore flags
            ])
            
            # New group
            group = dict()
            
            # Read headers
            for line in file:
                # Check if we found the correct header
                match = re.match(header[0][0], line)
                if not match:
                    raise IOError('Could not parse group header line: '+line)
                
                # Save matched expression?
                if header[0][1]:
                    group[header[0][1]] = header[0][2](match)
                    
                # Header found -> goto next header
                header.popleft()
                        
                # All headers -> stop
                if not header:
                    break
                
            # Split header information
            group['size'] = group['header'][0]
            group['material'] = group['header'][1]
            del group['header']
                
            # Add cell list
            group['cells'] = []
                
            # Read cells
            for line in file:
                # End of section?
                if line == ENDSECTION:
                    break
                
                cells = line.split()
                if len(cells) > 10:
                    raise IOError('Too many cells in one line in group section: '+line)
                
                group['cells'].extend(map(lambda x: int(x)-1, cells))
                
            # Check group size
            if group['size'] != len(group['cells']):
                raise IOError('Size of group '+group['name']+' not correct: size in header = '
                    +str(group['size'])+' != actual size = '+str(len(group['cells'])))
                
            # Add group
            self.__groups.append(group)
            
        else:
            # Add dummy group
            self.__groups.append(None)
            file.skipSection()
        
    def __readBoundary(self, file, save):
        if save:
            # Very strict at the moment, we might need to relax this
            header = re.compile('^(\d+) +1 +(\d+) +0 +6$')
            value = re.compile('^(\d+) +6 +(\d+)')
            
            # New boundary
            boundary = dict()
            
            # Read header
            for line in file:
                # Check for correct header
                match = re.match(header, line)
                if not match:
                    raise IOError('Could not parse boundary header line: '+line)
                
                # Get values
                boundary['name'] = match.group(1)
                boundary['size'] = int(match.group(2))
                
                # We only have one header line
                break
            
            boundary['sides'] = []
            
            # Read values
            for line in file:
                # Check for boundary
                match = re.match(value, line)
                if match:
                    boundary['sides'].append([int(match.group(1))-1, int(match.group(2))])
                    continue
                
                # End of section?
                if line == ENDSECTION:
                    break;
                
                raise IOError('Error in cell section in line: '+line)
            
            # Check size
            if boundary['size'] != len(boundary['sides']):
                raise IOError('Size of boundary '+boundary['name']+' not correct: size in header = '
                    +str(boundary['size'])+' != actual size = '+str(len(boundary['sides'])))
                
            # Add new boundary type
            self.__boundaries.append(boundary)
            
        else:
            # Add dummy boundary
            self.__boundaries.append(None)
            file.skipSection()
        
    def coords(self):
        return self.__coords
        
    def elements(self):
        return self.__cells
    
    def groups(self):
        return self.__groups
    
    def boundaries(self):
        return self.__boundaries
    
class GambitWriter:
    """Writes a GAMBIT mesh"""
    
    def __init__(self, mesh, file):
        # Open file if a string is specified
        if isinstance(file, basestring):
            file = open(file, 'w')
            
        minElementIdSize = int(math.floor(math.log10(len(mesh.elements())))) + 1
        minVertexIdSize = int(math.floor(math.log10(len(mesh.coords())))) + 1
            
        # Header
        print >> file, '%20s' % ('CONTROL INFO'), mesh.version
        print >> file, '** GAMBIT NEUTRAL FILE'
        print >> file, mesh.name
        print >> file, 'PROGRAM:  %20s     VERSION:  %s' % ('Gambit',   mesh.version)
        print >> file, ' ' + mesh.date
        print >> file, ' %9s %9s %9s %9s %9s %9s' % ('NUMNP', 'NELEM', 'NGRPS', 'NBSETS', 'NDFCD', 'NDFVL')
        print >> file, ' %9d %9d %9d %9d %9d %9d' % (len(mesh.coords()), len(mesh.elements()),
            len(mesh.groups()), len(mesh.boundaries()), 3, 3)
        print >> file, ENDSECTION
        
        # Coords
        print >> file, '%20s' % ('NODAL COORDINATES'), mesh.version
        for i, coord in enumerate(mesh.coords()):
            print >> file, '%10d %19.10e %19.10e %19.10e' % (i+1, coord[0], coord[1], coord[2])
        print >> file, ENDSECTION
        
        # Elements
        elementIdSize = str(max(7, minElementIdSize))
        vertexIdSize = str(max(7, minVertexIdSize))
        print >> file, '%20s' % ('ELEMENTS/CELLS'), mesh.version
        for i, element in enumerate(mesh.elements()):
            print >> file, (' %'+elementIdSize+'d %2d %2d' \
                + '  %'+vertexIdSize+'d %'+vertexIdSize+'d %'+vertexIdSize+'d %'+vertexIdSize+'d') \
                % ((i+1, 6, 4) + tuple([x+1 for x in element]))
        print >> file, ENDSECTION
        
        # Groups
        elementIdSize = str(max(7, minElementIdSize+1))
        for i, group in enumerate(mesh.groups()):
            print >> file, '%20s' % ('ELEMENT GROUP'), mesh.version
            print >> file, 'GROUP: %10d ELEMENTS: %10d MATERIAL: %10d NFLAGS: %10d' % (i+1, group['size'], group['material'], 1)
            print >> file, '%32s' % (group['name'])
            print >> file, '%8d' % (0)
            for i, cell in enumerate(group['cells']):
                # Whitespace at beginning of each line
                if i % 10 == 0:
                    print >> file, '',
                print >> file, ('%'+elementIdSize+'d') % (cell+1),
                # Newline after 10 elements or when finished
                if (i+1) % 10 == 0 or i+1 == len(group['cells']):
                    print >> file
            print >> file, ENDSECTION
            
        # Boundary conditions
        for boundary in mesh.boundaries():
            print >> file, '%20s' % ('BOUNDARY CONDITIONS'), mesh.version
            print >> file, '%32s%8d%8d%8d%8d' % (boundary['name'], 1, boundary['size'], 0, 6)
            for side in boundary['sides']:
                print >> file, '%10d %4d %4d' % (side[0]+1, 6, side[1])
            print >> file, ENDSECTION
            
        file.close()
