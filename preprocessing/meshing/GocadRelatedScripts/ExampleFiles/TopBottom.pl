#Author Thomas Ulrich
#Gocad Example for creating 2 plines to delimitate top and bottom faces of the domain
#more information about GOCAD formats can be found here
#http://paulbourke.net/dataformats/gocad/gocad.pdf
GOCAD Pline 1
HEADER {
name:bottom
}
ILINE
VRTX 1 60e3 60e3 -60e3
VRTX 2 -60e3 60e3 -60e3
VRTX 3 -60e3 -60e3 -60e3
VRTX 4 60e3 -60e3 -60e3
SEG 1 2 
SEG 2 3 
SEG 3 4
SEG 4 1
END
GOCAD Pline 1
HEADER {
name:top
}
ILINE
VRTX 1 60e3 60e3 -0e3
VRTX 2 -60e3 60e3 -0e3
VRTX 3 -60e3 -60e3 -0e3
VRTX 4 60e3 -60e3 -0e3
SEG 1 2 
SEG 2 3 
SEG 3 4
SEG 4 1
END
