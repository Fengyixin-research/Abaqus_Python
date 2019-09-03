#coding:utf-8
### Output nodes of the original model
from abaqusConstants import *
from odbAccess import *
odb = openOdb(path = 'Job-MultiPlateContact.odb')
assembly = odb.rootAssembly
# First step of the simulation
step1 = odb.steps['Step-ThermalExpansion']
lastFrame = step1.frames[-1]
# Number of deformable elements in X- and Y-directions
nX = 5
nY = 5
# Obtain the entire node set
for i in range(nX+1):
    for j in range(nY+1):
        nodeName = 'NODE-PLATE-'+str(i)+'-'+str(j)
        nodes = assembly.nodeSets[nodeName]
        displacement = lastFrame.fieldOutputs['U']
        nodesDisplacement = displacement.getSubset(region = nodes)
        # Output node information
        for v in nodesDisplacement.values:
            print 'Node Number: ', v.nodeLabel
            print 'Magnitude: ', v.magnitude
            print ' X-displacement:', v.data[0], ' Y-displacement:', v.data[1]
odb.close()

