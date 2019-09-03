#coding:utf-8
### Output node of the original model
from abaqusConstants import *
from odbAccess import *
odb = openOdb(path = 'Job-MultiPlateContact.odb')
assembly = odb.rootAssembly
# First step of the simulation
step1 = odb.steps['Step-ThermalExpansion']
lastFrame = step1.frames[-1]
# Obtain the node set containing the origin point
nodes = assembly.nodeSets['SET-PLATE-ORIGIN-FIX']
displacement = lastFrame.fieldOutputs['U']
nodesDisplacement = displacement.getSubset(region = nodes)
# Output node information
for v in nodesDisplacement.values:
    print 'Location: ', v.position
    print 'Type: ', v.type
    print 'Node Number: ', v.nodeLabel, ' X-displacement: ', v.data[0], ' Y-displacement: ', v.data[1]
    print 'X-displacement: ', v.data[0]
    print 'Y-displacement: ', v.data[1]
    print 'Magnitude: ', v.magnitude
odb.close()


