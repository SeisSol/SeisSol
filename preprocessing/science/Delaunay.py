##
# @file
# This file is part of SeisSol.
#
# @section LICENSE
# Copyright (c) SeisSol Group
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

import numpy as np
from scipy.spatial import Delaunay
import os
import timeit


class CheckStationsMesh(object):

    def __init__(
        self, filename_mesh='/home/msimon/svn/repos/verce/All/JRA/JRA1/python/test_ressources/inputfiles/eucrust_small_new.neu',
            station_coordinates=str()):
        """
        Checks if stations given by station_coordinates are inside the mesh.
        If not, the function "move_points_to_mesh" will move the points
        into the mesh and return a numpy ndarray of the coordinates of
        all the stations that were successfully placed into the mesh

        :type filename: class:'~str'
        :param: Fullpath, or relative path of the mesh
        :type station_coordinates: class: '~str'
        :param: String containing station coordinats as in Input.par
        e.g.:
          '3.9463538e+06   1.9802029e+06   4.5486094e+06 \n
           1.9542831e+06   9.9590736e+05   5.9475798e+06 \n
           4.1565145e+06   1.4021579e+06   4.5762539e+06
           .
           .
           .'

        """

        self.filename_mesh = filename_mesh
        try:
            self.station_coordinates = self.__read_station_coordinates(
                station_coordinates)
        except:
            print 'cannot load station coordinates'
            raise Exception, ('cannot load station coordinates %s' %
                              station_coordinates)
        try:
            self.vertices, self.elements, self.boundary_elements = \
                self.__read_mesh(filename_mesh=self.filename_mesh)
        except:
            print 'cannot load mesh'
            raise Exception, ('cannot load mesh %s' % self.filename_mesh)
        try:
            self.surface_volume = self.__construct_surface_volume(
                self.vertices, self.elements, self.boundary_elements,
                surface_volume_thickness=0.97)
        except:
            print 'cannot load construct surface volume'
            raise Exception('cannot construct surface volume')

    def __read_station_coordinates(self, station_coordinates):
        """
        Converts a String of station coordinates as given in Seissol-Input.par
        Input-file into a numpy array.

        :type station_coordinates: class:'str()'
        :param: Converts
        :rtype points: class:'numpy.ndarray()'
        :param: 3xN-dim array containing X,Y,Z coordinates of N stations
        """

        station_coordinates = station_coordinates.split('\n')
        points = []
        for station in station_coordinates:
            if station in '':
                break
            points.append([float(coord) for coord in station.split()])
        points = np.asanyarray(points)
        return points

    def __read_mesh(self, filename_mesh='/home/msimon/svn/repos/verce/All/JRA/JRA1/python/test_ressources/inputfiles/eucrust_small_new.neu'):
        """
        From mesh-file read vertices, elements and boundary elements

        :type: filename_mesh: 'str'
        :param: filename of the mesh. Full path or relative path
        :rtype: vertices: '~numpy.ndarray'
        :param: 4xN-dim array (N=Number of vertices). 1st-dim=Index;
        2nd-dim= x-coordinate; 3rd-dim = y-coordinate;
        4th-dim = z-coordinate
        :rtype: elements: '~numpy.dnarray'
        :param: 7xN-dim array (N=Number of elements). 1st-dim=Index;
        2nd-dim= irrelevant; 3rd-dim = irrelevant;
        4th-7th-dim = Indices of vertices constructing the element
        :rtype: boundary_elements: 'dict(boundary_condition: numpy.array)'
        :param: Dictionary Refering to the Indices of Elements adhering \
        to boundary conditions. Keys are Boundary conditions type
        identifiers
        Identifier of boundary element type:
        101 = free surface 	    boundary
        102 = non-conforming 	boundary
        103 = dynamic rupture	boundary
        104 = inflow 		    boundary
        105 = absorbing 	 	boundary
        106 = periodic 		    boundary
        """

        filename = filename_mesh
        if filename[-4:] not in '.neu':
            print 'filename not in *.neu format!!'
        filename = filename[:-4]
        t = timeit.time.time()
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #%        Read Gambit Data
        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f = open(filename + '.neu')
        print '\n'
        print '--------------------------------------------------------',\
            '---------------------------'
        print ' Reading data from: %s' % filename + '.neu'
        for i in range(6):
            junk = f.readline()
        tmp = f.readline()
        tmp = [int(s) for s in tmp.split() if s.isdigit()]
        NUMNP = tmp[0]  # %Number of Vertices
        NELEM = tmp[1]  # %Number of Elements
        NGRPS = tmp[2]  # %Number of Element Groups
        NBSETS = tmp[3]  # %Number of Boundary Condition Sets
        NDFCD = tmp[4]  # %Number of Dimensions (2 or 3)
        if 3 == NDFCD:

            # Vertices
            f.readline()
            f.readline()
            vertices = np.fromfile(file=f, count=NUMNP * 4, sep=" ")
            vertices = vertices.reshape(NUMNP, 4)
            vertices = vertices.transpose()

            # Elements
            f.readline()
            f.readline()
            position = f.tell()
            first_line = f.readline()
            first_line = [int(s) for s in first_line.split() if s.isdigit()]
            n = 7
            f.seek(position, 0)
            elements = np.fromfile(file=f, count=NELEM * n, sep=" ",
                                   dtype=int).reshape(NELEM, n).transpose()

            # loop over element groups. irrelevant for this routine, but i\
            # left it in in case someone wants to work with it
            for group in range(NGRPS):
                f.readline()
                f.readline()
                first_line = f.readline()
                first_line = [int(s) for s in first_line.split() if
                              s.isdigit()]
                nelem = first_line[1]
                f.readline()
                f.readline()
                tmp = np.fromfile(file=f, count=nelem, sep=" ")
                position = f.tell()

            # Boundary Elements
            f.seek(position, 0)
            boundary_elements = {}
            for boundaries in range(NBSETS):
                f.readline()
                f.readline()
                first_line = f.readline()
                first_line = [int(s) for s in first_line.split() if
                              s.isdigit()]
                nelem = first_line[2]
                boundary_elements[str(
                    first_line[0])] = np.fromfile(
                        file=f, count=nelem * 3, sep=" ", dtype=int).reshape(3,
                                                                             nelem)[0]
                position = f.tell()
        else:
            print 'only 3d tetrahedral meshes'
        self.vertices = vertices
        self.elements = elements
        self.boundary_elements = boundary_elements
        return vertices, elements, boundary_elements

    def __construct_surface_volume(self, vertices, elements,
                                   boundary_elements, surface_volume_thickness=0.05):
        """
        From Free surface boundary construct convex hull of
        surface volume

        :type vertices: class:'~numpy.ndarray'
        :param: 4xN-dim array (N=Number of vertices). 1st-dim=Index;
        2nd-dim= x-coordinate; 3rd-dim = y-coordinate;
        4th-dim = z-coordinate
        :type elements: class:'~numpy.dnarray'
        :param: 7xN-dim array (N=Number of elements). 1st-dim=Index;
        2nd-dim= irrelevant; 3rd-dim = irrelevant;
        4th-7th-dim = Indices of vertices constructing the element
        :type boundary_elements: class:'dict(boundary_condition: numpy.array)'
        :param: Dictionary Refering to the Indices of Elements adhering \
        to boundary conditions.
        :type surface_volume_thickness: class:'~float'
        :param: Given in percent. Defines the volume for which to
        calculate the convex hull below the surface.
        :rtype surface_volume: class:'~numpy.ndarray'
        :param: Collection of point coordinates of surface element-vertices
        and the same points shifted towards coordinate system origin by
        some percent defined in surface_volume_thickness
        """

        # select elements for boundary condition '101' = free surface
        boundary_vertices = elements[3:].transpose().take(
            boundary_elements['101']-1, axis=0)
        # extract coordinates of vertices
        boundary_points_coordinates = vertices[1:].transpose().take(
            boundary_vertices.flatten()-1, axis=0)
        # create a surface volume by adding the same points, just moved
        # slightly towards the coordinate system origin
        lower_bound_boundary_points = boundary_points_coordinates *\
            (1.0 - surface_volume_thickness)
        surface_volume = np.concatenate((boundary_points_coordinates,
                                         lower_bound_boundary_points), axis=0)
        return surface_volume

    def __in_hull(self, points, hull):
        """
        hull is Delaunay of scipy.spatial

        :type points: class:'~numpy.ndarray'
        :params: 3xN-dim numpy array containing the X,Y,Z coordinates of the
        N stations.
        :type hull: class:'~scipy.spatial.Delaunay'
        :params: convex hull around mesh surface
        :rtype: 1d-Boolean array stating if Point is inside the convex
        hull or not
        """

        return hull.find_simplex(points) >= 0

    def move_points_to_mesh(self, steplength=1.01):
        """
        Checks if the stations are located within the mesh.

        The stations outside the mesh will be moved incrementally towards
        the coordinate system origin, until they are in the mesh again.
        The initial position of each station is moved away from the origin
        by a factor of "steplength"
        Then moving back towards the origin, is done in 10ths of steplength

        If some coordinates could not be moved into the mesh, they will
        be ignored.

        :type steplength: class:'~float'
        :param: sets the inital offset of the station coordinates away
        from the origin
        :rtype points_in_mesh: class:'~numpy.ndarray'
        :param: 3xN-dim numpy array containing all coordinates that could
         successfully be moved into the mesh.
        """
        
      
        
        try:
            self.hull = Delaunay(self.surface_volume)
        except:
            raise Exception, "Cant create convex hull for mesh"

        # set inital values:
        finalstation_coordinates = np.copy(self.station_coordinates *
                                           steplength)
        inhull = self.__in_hull(finalstation_coordinates, self.hull)
        min_radius_finalstation_coordinates = 1.0
        min_radius_surface_volume = 0.0

        # move stations into the convex hull of the surface volume. stop,
        # in case their distance to the origin is smaller then the the
        # bottom layer of the surface volume
        while min_radius_finalstation_coordinates > \
            min_radius_surface_volume and \
                sum(inhull) != len(self.station_coordinates):
            inhull = self.__in_hull(finalstation_coordinates, self.hull)
            inv_in3d = np.column_stack((inhull, inhull, inhull))
            finalstation_coordinates = np.where(inv_in3d,
                                                finalstation_coordinates, finalstation_coordinates *
                                               (1 - steplength / 10))
            min_radius_finalstation_coordinates = min([np.sqrt(x * x + y * y + z * z)
                                                       for x, y, z
                                                       in finalstation_coordinates])
            min_radius_surface_volume = min([np.sqrt(x * x + y * y + z * z) for
                                             x, y, z
                                             in self.surface_volume])

        # only keep the points that could be moved into the mesh
        points_in_mesh = np.compress(inv_in3d.flatten(),
                                     finalstation_coordinates.
                                     flatten(
                                     )).reshape(
                                         len(np.compress(inv_in3d.flatten(),
                                                         finalstation_coordinates.flatten())) / 3, 3)

        return points_in_mesh
        

if __name__ == '__main__':
    
    f=open('/home/msimon/svn/repos/verce/All/JRA/JRA1/python/'+\
            'test_ressources/inputfiles/input.par')
    strin=f.read()
    station_coordinates='\n'.join(strin.split('\n')[88:-21])
    #csm = CheckStationsMesh(filename_mesh='/home/msimon/svn/repos/verce/All/JRA/JRA1/python/\
    #test_ressources/inputfiles/eucrust_small_new.neu', \
    csm = CheckStationsMesh(filename_mesh='/home/msimon/git/'+\
    'wfs_input_generator_msimon00/wfs_input_generator/tests/data/'+\
    'seissol_example/most_simple_tet.neu', \
    station_coordinates=station_coordinates)
    
    f.close()
    points_in_mesh = csm.move_points_to_mesh()
    
    print points_in_mesh
