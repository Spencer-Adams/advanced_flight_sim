import numpy as np
from stl import mesh
import math as m
PRECISION = 1.0e-5 # points within this distance will be treated as non-unique
stlfile = 'F16_coarse.stl' # input file
vtkfile = 'F16_coarse.vtk' # output file
myMesh = mesh.Mesh.from_file(stlfile)
meshpoints_orig = myMesh.vectors
nfaces = len(myMesh.vectors)
# Get unique lines
templines = np.zeros((3*nfaces,2,3)) # initialize first two lines
print("Removing duplicate lines for ",nfaces," faces.")
num_unique = 0
for f in range(nfaces): # loop over all the faces
    print("Checking face ",f," of ",nfaces," faces")
    for i in range(3): # loop over all the lines on each face
        if(i+1 < 3): # if looking at last vertex, use vertex 2 and 0 for comparison
            ip1 = i+1
        else:
            ip1 = 0
        pa = myMesh.vectors[f,i,:]
        pb = myMesh.vectors[f,ip1,:]
        found = False
        for j in range(num_unique): # loop through the unique lines in lines3D
            p1 = templines[j,0,:]
            p2 = templines[j,1,:]
            dx1 = pa[0] - p1[0]
            dy1 = pa[1] - p1[1]
            dz1 = pa[2] - p1[2]
            dx2 = pb[0] - p2[0]
            dy2 = pb[1] - p2[1]
            dz2 = pb[2] - p2[2]
            if(abs(dx1)+abs(dy1)+abs(dz1)+abs(dx2)+abs(dy2)+abs(dz2) < PRECISION):
                found = True #lines match
                break
            else: # check to see if other direction matches
                dx1 = pa[0] - p2[0]
                dy1 = pa[1] - p2[1]
                dz1 = pa[2] - p2[2]
                dx2 = pb[0] - p1[0]
                dy2 = pb[1] - p1[1]
                dz2 = pb[2] - p1[2]
            if(abs(dx1)+abs(dy1)+abs(dz1)+abs(dx2)+abs(dy2)+abs(dz2) < PRECISION):
                found = True #lines match
                break
        if(not found): #lines don't match, so add lines to unique list
            templines[num_unique,0,:] = pa
            templines[num_unique,1,:] = pb
            num_unique += 1

nlines = num_unique
print("Number of unique lines = ",nlines)
lines = np.zeros((nlines,2,3))
lines = templines[0:nlines,:,:]
# Get unique points
points = np.zeros((1,3))
points[0,:] = lines[0,0,:]
lineid = np.zeros((nlines,2), dtype=int)
print("Removing duplicate points for ",nlines," lines.")
num_unique = 1

for f in range(nlines): # loop over all the lines
    print("Checking line ",f," of ",nlines," lines")
    for i in range(2): # loop over all the points on each line
        pa = lines[f,i,:]
        found = False
        for j in range(num_unique): # loop through the unique points
            p1 = points[j,:]
            dx = pa[0] - p1[0]
            dy = pa[1] - p1[1]
            dz = pa[2] - p1[2]
            if(abs(dx) + abs(dy) + abs(dz) < PRECISION): #points match
                found = True
                lineid[f,i] = j
                break
        if(not found): #points don't match, so add point to unique list
            points = np.append(points,[pa],axis=0)
            lineid[f,i] = num_unique
            num_unique += 1
npoints = len(points)
# Create VTK file and write only unique points
with open(vtkfile, 'w') as writer:
    writer.write('# vtk DataFile Version 3.0\n')
    writer.write(f'{vtkfile}\n')
    writer.write('ASCII\n')
    writer.write('DATASET POLYDATA\n')
    writer.write(f'POINTS {npoints} float\n')
    for i in range(npoints):
        writer.write(f'{points[i,0]:.6e} {points[i,1]:.6e} {points[i,2]:.6e}\n')
    writer.write(f'LINES {nlines} {3*nlines}\n')   
    for i in range(nlines):
        writer.write(f'2 {lineid[i,0]} {lineid[i,1]}\n')

print('')
print("Finished writing ",vtkfile)
print('')
print("File Statistics:")
print("--> Original number of points = ",nfaces*3)
print("--> Original number of lines = ",nfaces*3)
print("--> VTK number of unique points = ",npoints)
print("--> VTK number of unique lines = ",nlines)

    



  
