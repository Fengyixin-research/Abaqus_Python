#coding:utf-8
from abaqus import *
import testUtils
testUtils.setBackwardCompatibility()
from abaqusConstants import *
import __main__
import section, regionToolset, connectorBehavior, random
import displayGroupMdbToolset as dgm
import part, material, assembly, step, interaction, load, mesh, optimization, job, sketch, visualization, xyPlot
import displayGroupOdbToolset as dgo
from odbAccess import *
from deap import base
from deap import creator
from deap import tools
import csv, os
## Import coordinates information from files
scaling = 1000.0
dataUpper, dataLower = [], []
# Import Upper3
with open('Upper3.csv','r') as csvfile:
    csv_reader = csv.reader(csvfile) 
    for row in csv_reader:
        dataUpper.append(row)
# Import Lower3
with open('Lower3.csv','r') as csvfile:
    csv_reader = csv.reader(csvfile) 
    for row in csv_reader:
        dataLower.append(row)
UpperXY, LowerXY = [], []
for item in dataUpper:
    dataxy=[scaling*float(xy) for xy in item]
    UpperXY.append(dataxy)
for item in dataLower:
    dataxy=[scaling*float(xy) for xy in item]
    LowerXY.append(dataxy)
# Point left--->right
LowerXY=LowerXY[::-1]

#############################################
## Number of elements in each dimension
nX, nY = len(UpperXY)-1, len(UpperXY[0])/2-1
nXL, nYL = len(LowerXY)-1, len(LowerXY[0])/2-1
lengthMesh = 0.005*scaling
# Set temperature for initial and final state
initialT, finalT = 0.0, 100.0
numDomainsCPUs = 4
sizePopulation, numIter, nsize = 30, 25, 5
expandActive, expandInactive = 0.0005, 0
offset = 1.5
lower, upper = 0, 1
## Upper and Lower origin point coordinates
UpperOrigin=[UpperXY[0][0], UpperXY[0][1]]
LowerOrigin=[LowerXY[0][0], LowerXY[0][1]]
#############################################

############### Calculation of objective function ############### 
# 3 parameters: individual object, iteration number, model number
def evalOneMax(individual, iter, nmodel):
    ## Create a new model named 'Model'
    modelName = 'Model-Iter'+str(iter)+'-Num'+str(nmodel)
    myModel = mdb.Model(name= modelName)
    ## Create Upper constraint sketches, starting from the top-left corner
    for i in range(nX):
        for j in range(nY):
            s = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=scaling)
            g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
            s.setPrimaryObject(option=STANDALONE)
            # Points in the data array is sorted according to decreasing Y-coordinates
            s.Line(point1=(UpperXY[i][2*j+0], UpperXY[i][2*j+1]),     point2=(UpperXY[i+1][2*j+0], UpperXY[i+1][2*j+1]))
            s.Line(point1=(UpperXY[i+1][2*j+0], UpperXY[i+1][2*j+1]), point2=(UpperXY[i+1][2*j+2], UpperXY[i+1][2*j+3]))
            s.Line(point1=(UpperXY[i+1][2*j+2], UpperXY[i+1][2*j+3]), point2=(UpperXY[i][2*j+2], UpperXY[i][2*j+3]))
            s.Line(point1=(UpperXY[i][2*j+2], UpperXY[i][2*j+3]),     point2=(UpperXY[i][2*j+0], UpperXY[i][2*j+1]))
            # Create part named 'Part-Upper-i-j'
            partName = 'Part-Upper-'+str(i)+'-'+str(j)
            p = mdb.models[modelName].Part(name=partName, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
            p = mdb.models[modelName].parts[partName]
            p.BaseShell(sketch=s)
            s.unsetPrimaryObject()
            del mdb.models[modelName].sketches['__profile__']
    ## Create Lower constraint sketches, starting from the bottom-left corner
    for i in range(nXL):
        for j in range(nYL):
            s = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=scaling)
            g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
            s.setPrimaryObject(option=STANDALONE)
            # Points in the data array is sorted according to increasing Y-coordinates
            s.Line(point1=(LowerXY[i][2*j+0], LowerXY[i][2*j+1]),     point2=(LowerXY[i+1][2*j+0], LowerXY[i+1][2*j+1]))
            s.Line(point1=(LowerXY[i+1][2*j+0], LowerXY[i+1][2*j+1]), point2=(LowerXY[i+1][2*j+2], LowerXY[i+1][2*j+3]))
            s.Line(point1=(LowerXY[i+1][2*j+2], LowerXY[i+1][2*j+3]), point2=(LowerXY[i][2*j+2], LowerXY[i][2*j+3]))
            s.Line(point1=(LowerXY[i][2*j+2], LowerXY[i][2*j+3]),     point2=(LowerXY[i][2*j+0], LowerXY[i][2*j+1]))
            # Create part named 'Part-Lower-i-j'
            partName = 'Part-Lower-'+str(i)+'-'+str(j)
            p = mdb.models[modelName].Part(name=partName, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
            p = mdb.models[modelName].parts[partName]
            p.BaseShell(sketch=s)
            s.unsetPrimaryObject()
            del mdb.models[modelName].sketches['__profile__']
    ## Create tip of the wing
    s = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=scaling)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    for j in range(nY):
        s.Line(point1=(UpperXY[nX][2*j+0], UpperXY[nX][2*j+1]),   point2=(UpperXY[nX][2*(j+1)+0], UpperXY[nX][2*(j+1)+1]))
        s.Line(point1=(LowerXY[nXL][2*j+0], LowerXY[nXL][2*j+1]), point2=(LowerXY[nXL][2*(j+1)+0], LowerXY[nXL][2*(j+1)+1]))
    s.Line(point1=(UpperXY[nX][0], UpperXY[nX][1]),   point2=(scaling*1.0, 0.0))
    s.Line(point1=(LowerXY[nXL][0], LowerXY[nXL][1]), point2=(scaling*1.0, 0.0))
    partName = 'Part-Tip'
    p = mdb.models[modelName].Part(name=partName, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    p = mdb.models[modelName].parts[partName]
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()
    del mdb.models[modelName].sketches['__profile__']

    ## Create materials named 'Material-Active' and 'Material-Inactive'
    activeName = 'Material-Active'
    inactiveName = 'Material-Inactive' 
    mdb.models[modelName].Material(name=activeName)
    mdb.models[modelName].materials[activeName].Elastic(table=((210000.0, 0.3), ))
    mdb.models[modelName].materials[activeName].Expansion(table=((expandActive, ), ))
    mdb.models[modelName].Material(name=inactiveName)
    mdb.models[modelName].materials[inactiveName].Elastic(table=((210000.0, 0.3), ))
    mdb.models[modelName].materials[inactiveName].Expansion(table=((expandInactive, ), ))
    ## Create HomogeneousSolidSection objects and assign material property to them 
    activeSection = 'Section-Active'
    inactiveSection = 'Section-Inactive'
    mdb.models[modelName].HomogeneousSolidSection(name=activeSection, material=activeName, thickness=None)
    mdb.models[modelName].HomogeneousSolidSection(name=inactiveSection, material=inactiveName, thickness=None)
    ## Assign material sections to Upper part elements
    for i in range(nX):
        for j in range(nY):
            # Create partName named 'Part-Upper-i-j'
            partName = 'Part-Upper-'+str(i)+'-'+str(j)
            p = mdb.models[modelName].parts[partName]
            f = p.faces
            faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
            region = regionToolset.Region(faces=faces)
            if individual[j*nX+i] == 0:
                p.SectionAssignment(region=region, sectionName=inactiveSection, 
                    offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
            elif individual[j*nX+i] == 1:
                p.SectionAssignment(region=region, sectionName=activeSection, 
                    offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
    ## Assign material sections to Lower part elements
    for i in range(nXL):
        for j in range(nYL):
            # Create partName named 'Part-Lower-i-j'
            partName = 'Part-Lower-'+str(i)+'-'+str(j)
            p = mdb.models[modelName].parts[partName]
            f = p.faces
            faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
            region = regionToolset.Region(faces=faces)
            if individual[j*nXL+i+nX*nY] == 0:
                p.SectionAssignment(region=region, sectionName=inactiveSection, 
                    offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
            elif individual[j*nXL+i+nX*nY] == 1:
                p.SectionAssignment(region=region, sectionName=activeSection, 
                    offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    partName='Part-Tip'
    p = mdb.models[modelName].parts[partName]
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    region = regionToolset.Region(faces=faces)
    p.SectionAssignment(region=region, sectionName=inactiveSection, offset=0.0, offsetType=MIDDLE_SURFACE, 
        offsetField='', thicknessAssignment=FROM_SECTION)

    ## Create assembly consisting of all plate instances
    a = mdb.models[modelName].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    for i in range(nX):
        for j in range(nY):
            partName = 'Part-Upper-'+str(i)+'-'+str(j)
            instanceName = partName+'-1'
            p = mdb.models[modelName].parts[partName]
            # Create indenpendent instances by setting 'dependent' to OFF
            a.Instance(name=instanceName, part=p, dependent=OFF)
    for i in range(nXL):
        for j in range(nYL):
            partName = 'Part-Lower-'+str(i)+'-'+str(j)
            instanceName = partName+'-1'
            p = mdb.models[modelName].parts[partName]
            a.Instance(name=instanceName, part=p, dependent=OFF)
    partName='Part-Tip'
    instanceName = partName+'-1'
    p = mdb.models[modelName].parts[partName]
    a.Instance(name=instanceName, part=p, dependent=OFF)

    ## Generate mesh for the assembly
    partInstances = []
    for i in range(nX):
        for j in range(nY):
            instanceName = 'Part-Upper-'+str(i)+'-'+str(j)+'-1'
            partInstances.append(a.instances[instanceName])
    for i in range(nXL):
        for j in range(nYL):
            instanceName = 'Part-Lower-'+str(i)+'-'+str(j)+'-1'
            partInstances.append(a.instances[instanceName])
    instanceName = 'Part-Tip-1'
    partInstances.append(a.instances[instanceName])
    partInstances=tuple(partInstances)
    a.seedPartInstance(regions=partInstances, size=lengthMesh, deviationFactor=0.1, minSizeFactor=0.1)
    a.generateMesh(regions=partInstances)
    ## Create a new step named 'Step-ThermalExpansion' after initial step
    mdb.models[modelName].StaticStep(name='Step-ThermalExpansion', previous='Initial')

    ################## Define interaction relations ##################
    ## Create Upper interaction edges and define interaction relations 
    for i in range(nX):
        for j in range(nY):
            # Create instanceName named 'Part-Upper-i-j-1'
            instanceName = 'Part-Upper-'+str(i)+'-'+str(j)+'-1'
            # Lower:3, Right:2, Upper:1, Left:4
            s1 = a.instances[instanceName].edges
            # Edges are located in clockwise direction
            edgeCenter1 = (0.5*(UpperXY[i][2*j+0]+UpperXY[i+1][2*j+0]),   0.5*(UpperXY[i][2*j+1]+UpperXY[i+1][2*j+1]), 0)
            edgeCenter2 = (0.5*(UpperXY[i+1][2*j+0]+UpperXY[i+1][2*j+2]), 0.5*(UpperXY[i+1][2*j+1]+UpperXY[i+1][2*j+3]), 0)
            edgeCenter3 = (0.5*(UpperXY[i+1][2*j+2]+UpperXY[i][2*j+2]),   0.5*(UpperXY[i+1][2*j+3]+UpperXY[i][2*j+3]), 0)
            edgeCenter4 = (0.5*(UpperXY[i][2*j+2]+UpperXY[i][2*j+0]),     0.5*(UpperXY[i][2*j+3]+UpperXY[i][2*j+1]), 0)
            edges1 = s1.findAt((edgeCenter1, ))
            edges2 = s1.findAt((edgeCenter2, ))
            edges3 = s1.findAt((edgeCenter3, ))
            edges4 = s1.findAt((edgeCenter4, ))
            edges1 = (edges1, ); edges2 = (edges2, ); edges3 = (edges3, ); edges4 = (edges4, )
            # Create surfaces for all 4 edges of current element
            a.Surface(side1Edges=edges1, name='Surf-Upper-'+str(i)+'-'+str(j)+'-1') 
            a.Surface(side1Edges=edges2, name='Surf-Upper-'+str(i)+'-'+str(j)+'-2') 
            a.Surface(side1Edges=edges3, name='Surf-Upper-'+str(i)+'-'+str(j)+'-3')  
            a.Surface(side1Edges=edges4, name='Surf-Upper-'+str(i)+'-'+str(j)+'-4')  
    # Define interaction relations between lower and upper surfaces
    for i in range(nX):
        for j in range(nY-1):
            surfName1 = 'Surf-Upper-'+str(i)+'-'+str(j)+'-3'
            surfName2 = 'Surf-Upper-'+str(i)+'-'+str(j+1)+'-1'
            region1=a.surfaces[surfName1]; region2=a.surfaces[surfName2]
            contactName = 'Constraint-Upper-'+str(i)+str(j)+str(i)+str(j+1)+'-31'
            mdb.models[modelName].Tie(name=contactName, master=region1, slave=region2, 
                positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
    # Define interaction relations between left and right surfaces
    for i in range(nX-1):
        for j in range(nY):
            surfName1 = 'Surf-Upper-'+str(i)+'-'+str(j)+'-2'
            surfName2 = 'Surf-Upper-'+str(i+1)+'-'+str(j)+'-4'
            region1=a.surfaces[surfName1]; region2=a.surfaces[surfName2]
            contactName = 'Constraint-Upper-'+str(i)+str(j)+str(i+1)+str(j)+'-24'
            mdb.models[modelName].Tie(name=contactName, master=region1, slave=region2, 
                positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

    ## Create Lower interaction edges and define interaction relations 
    for i in range(nXL):
        for j in range(nYL):
            # Create instanceName named 'Part-Plate-i-j-1'
            instanceName = 'Part-Lower-'+str(i)+'-'+str(j)+'-1'
            # Lower:1, Right:2, Upper:3, Left:4
            s1 = a.instances[instanceName].edges
            # Edges are located in clockwise direction
            edgeCenter1 = (0.5*(LowerXY[i][2*j+0]+LowerXY[i+1][2*j+0]),   0.5*(LowerXY[i][2*j+1]+LowerXY[i+1][2*j+1]), 0)
            edgeCenter2 = (0.5*(LowerXY[i+1][2*j+0]+LowerXY[i+1][2*j+2]), 0.5*(LowerXY[i+1][2*j+1]+LowerXY[i+1][2*j+3]), 0)
            edgeCenter3 = (0.5*(LowerXY[i+1][2*j+2]+LowerXY[i][2*j+2]),   0.5*(LowerXY[i+1][2*j+3]+LowerXY[i][2*j+3]), 0)
            edgeCenter4 = (0.5*(LowerXY[i][2*j+2]+LowerXY[i][2*j+0]),     0.5*(LowerXY[i][2*j+3]+LowerXY[i][2*j+1]), 0)
            edges1 = s1.findAt((edgeCenter1, ))
            edges2 = s1.findAt((edgeCenter2, ))
            edges3 = s1.findAt((edgeCenter3, ))
            edges4 = s1.findAt((edgeCenter4, ))
            edges1 = (edges1, ); edges2 = (edges2, ); edges3 = (edges3, ); edges4 = (edges4, )
            # Create surfaces for all 4 edges of current element
            a.Surface(side1Edges=edges1, name='Surf-Lower-'+str(i)+'-'+str(j)+'-1')  # Lower
            a.Surface(side1Edges=edges2, name='Surf-Lower-'+str(i)+'-'+str(j)+'-2')  # Left
            a.Surface(side1Edges=edges3, name='Surf-Lower-'+str(i)+'-'+str(j)+'-3')  # Upper
            a.Surface(side1Edges=edges4, name='Surf-Lower-'+str(i)+'-'+str(j)+'-4')  # Right
    # Define interaction relations between lower and upper surfaces
    for i in range(nXL):
        for j in range(nYL-1):
            surfName1 = 'Surf-Lower-'+str(i)+'-'+str(j)+'-3'
            surfName2 = 'Surf-Lower-'+str(i)+'-'+str(j+1)+'-1'
            region1=a.surfaces[surfName1]; region2=a.surfaces[surfName2]
            contactName = 'Constraint-Lower-'+str(i)+str(j)+str(i)+str(j+1)+'-31'
            mdb.models[modelName].Tie(name=contactName, master=region1, slave=region2, 
                positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
    # Define interaction relations between left and right surfaces
    for i in range(nXL-1):
        for j in range(nYL):
            surfName1 = 'Surf-Lower-'+str(i)+'-'+str(j)+'-2'
            surfName2 = 'Surf-Lower-'+str(i+1)+'-'+str(j)+'-4'
            region1=a.surfaces[surfName1]; region2=a.surfaces[surfName2]
            contactName = 'Constraint-Lower-'+str(i)+str(j)+str(i+1)+str(j)+'-24'
            mdb.models[modelName].Tie(name=contactName, master=region1, slave=region2, 
                positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
 
    ## Create interaction edges of the wing tip
    instanceName = 'Part-Tip-1'
    s1 = a.instances[instanceName].edges
    for j in range(nY):
        edgeCenter1 = (0.5*(UpperXY[nX][2*j+0]+UpperXY[nX][2*(j+1)+0]), 0.5*(UpperXY[nX][2*j+1]+UpperXY[nX][2*(j+1)+1]), 0)
        edgeCenter2 = (0.5*(LowerXY[nXL][2*j+0]+LowerXY[nXL][2*(j+1)+0]), 0.5*(LowerXY[nXL][2*j+1]+LowerXY[nXL][2*(j+1)+1]), 0)
        edges1 = s1.findAt((edgeCenter1, ))
        edges2 = s1.findAt((edgeCenter2, ))
        edges1 = (edges1, ); edges2 = (edges2, )
        a.Surface(side1Edges=edges1, name='Surf-Upper-Tip-'+str(j))
        a.Surface(side1Edges=edges2, name='Surf-Lower-Tip-'+str(j)) 
    # Define interaction relations
    for j in range(nY):
        surfName1 = 'Surf-Upper-'+str(nX-1)+'-'+str(j)+'-2'
        surfName2 = 'Surf-Upper-Tip-'+str(j)
        surfName3 = 'Surf-Lower-'+str(nXL-1)+'-'+str(j)+'-2'
        surfName4 = 'Surf-Lower-Tip-'+str(j)
        region1=a.surfaces[surfName1]
        region2=a.surfaces[surfName2]
        region3=a.surfaces[surfName3]
        region4=a.surfaces[surfName4]
        contactName1 = 'Constraint-Upper-Tip-'+str(j)
        contactName2 = 'Constraint-Lower-Tip-'+str(j)
        mdb.models[modelName].Tie(name=contactName1, master=region1, slave=region2, 
            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
        mdb.models[modelName].Tie(name=contactName2, master=region4, slave=region3, 
            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

    ## Create boundary conditions for the assembly
    a = mdb.models[modelName].rootAssembly
    for j in range(nY):
        instanceName = 'Part-Upper-0-'+str(j)+'-1'
        e1 = a.instances[instanceName].edges
        edgeCenter = (0.5*(UpperXY[0][2*j+0]+UpperXY[0][2*(j+1)+0]), 0.5*(UpperXY[0][2*j+1]+UpperXY[0][2*(j+1)+1]), 0.0)
        edges1 = e1.findAt((edgeCenter, )); edges1 = (edges1, )
        setName = 'Set-Upper-0-'+str(j)+'-Fix'
        a.Set(edges=edges1, name=setName)
        # Set boundary conditions
        region = a.sets[setName]
        bcName = 'BC-Upper-0-'+str(j)+'-Fix'
        mdb.models[modelName].DisplacementBC(name=bcName, createStepName='Initial', region=region, u1=SET, u2=UNSET, ur3=SET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    instanceName = 'Part-Upper-0-0-1'
    v1 = a.instances[instanceName].vertices
    vertexCenter = (UpperXY[0][0], UpperXY[0][1], 0.0)
    v = v1.findAt((vertexCenter, )); v = (v, )
    setName = 'Set-Upper-Origin-Fix'
    a.Set(vertices=v, name=setName)
    region = a.sets[setName]
    bcName = 'BC-Upper-Origin-Fix'
    # Origin fix in x-direction, y-direction, xy-rotation
    mdb.models[modelName].DisplacementBC(name=bcName, createStepName='Initial', region=region, u1=SET, u2=SET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

    ## Create boundary conditions for the assembly
    a = mdb.models[modelName].rootAssembly
    for j in range(nYL):
        instanceName = 'Part-Lower-0-'+str(j)+'-1'
        e1 = a.instances[instanceName].edges
        edgeCenter = (0.5*(LowerXY[0][2*j+0]+LowerXY[0][2*(j+1)+0]), 0.5*(LowerXY[0][2*j+1]+LowerXY[0][2*(j+1)+1]), 0.0)
        edges1 = e1.findAt((edgeCenter, )); edges1 = (edges1, )
        setName = 'Set-Lower-0-'+str(j)+'-Fix'
        a.Set(edges=edges1, name=setName)
        # Set boundary conditions
        region = a.sets[setName]
        bcName = 'BC-Lower-0-'+str(j)+'-Fix'
        mdb.models[modelName].DisplacementBC(name=bcName, createStepName='Initial', region=region, u1=SET, u2=UNSET, ur3=SET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
    instanceName = 'Part-Lower-0-0-1'
    v1 = a.instances[instanceName].vertices
    vertexCenter = (LowerXY[0][0], LowerXY[0][1], 0.0)
    v = v1.findAt((vertexCenter, )); v = (v, )
    setName = 'Set-Lower-Origin-Fix'
    a.Set(vertices=v, name=setName)
    region = a.sets[setName]
    bcName = 'BC-Lower-Origin-Fix'
    # Origin fix in x-direction, y-direction, xy-rotation
    mdb.models[modelName].DisplacementBC(name=bcName, createStepName='Initial', region=region, u1=SET, u2=SET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

    ## Define nodeSets that records the displacement information of the entire assembly
    ############# Upper #############
    # 1st: top-right vertices of all elements in the upper elements of the entire assembly
    a = mdb.models[modelName].rootAssembly
    for i in range(nX):
        nodeName = 'Node-Upper-'+str(i)+'-0'
        instanceName = 'Part-Upper-'+str(i)+'-0-1'
        v1 = a.instances[instanceName].vertices
        vertexCenter = (UpperXY[i+1][0], UpperXY[i+1][1], 0.0)
        v = v1.findAt((vertexCenter, )); v = (v, )
        a.Set(vertices=v, name=nodeName)
    # 2nd: leftmost vertex of all upper elements
    nodeName = 'Node-Upper-Origin'
    instanceName = 'Part-Upper-0-0-1'
    v1 = a.instances[instanceName].vertices
    vertexCenter = (UpperXY[0][0], UpperXY[0][1], 0.0)
    v = v1.findAt((vertexCenter, )); v = (v, )
    a.Set(vertices=v, name=nodeName)
    ############# Lower #############
    # 1st: bottom-right vertices of all elements in the lower elements of the entire assembly
    for i in range(nXL):
        nodeName = 'Node-Lower-'+str(i)+'-0'
        instanceName = 'Part-Lower-'+str(i)+'-0-1'
        v1 = a.instances[instanceName].vertices
        vertexCenter = (LowerXY[i+1][0], LowerXY[i+1][1], 0.0)
        v = v1.findAt((vertexCenter, )); v = (v, )
        a.Set(vertices=v, name=nodeName)
    # 2nd: leftmost vertex of all lower elements
    nodeName = 'Node-Lower-Origin'
    instanceName = 'Part-Lower-0-0-1'
    v1 = a.instances[instanceName].vertices
    vertexCenter = (LowerXY[0][0], LowerXY[0][1], 0.0)
    v = v1.findAt((vertexCenter, )); v = (v, )
    a.Set(vertices=v, name=nodeName)

    ## Define ambient temperature field
    a = mdb.models[modelName].rootAssembly
    allFaceSet = []
    for i in range(nX):
        for j in range(nY):
            instanceName = 'Part-Upper-'+str(i)+'-'+str(j)+'-1'
            f1 = a.instances[instanceName].faces
            faces = f1.getSequenceFromMask(mask=('[#1 ]', ), )
            allFaceSet.append(faces)
    for i in range(nXL):
        for j in range(nYL):
            instanceName = 'Part-Lower-'+str(i)+'-'+str(j)+'-1'
            f1 = a.instances[instanceName].faces
            faces = f1.getSequenceFromMask(mask=('[#1 ]', ), )
            allFaceSet.append(faces)
    instanceName = 'Part-Tip-1'
    f1 = a.instances[instanceName].faces
    faces = f1.getSequenceFromMask(mask=('[#1 ]', ), )
    allFaceSet.append(faces)
    allFaceSet =tuple(allFaceSet)
    setName = 'Set-UpperLower-AllFace'
    a.Set(faces=allFaceSet, name=setName)
    region = a.sets[setName]
    mdb.models[modelName].Temperature(name='Predefined Field-Temperature', createStepName='Initial', region=region, 
        distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(initialT, ))
    mdb.models[modelName].predefinedFields['Predefined Field-Temperature'].setValuesInStep(
        stepName='Step-ThermalExpansion', magnitudes=(finalT, ))

    ## Create job and submit. Set parallel computing
    jobName = 'Job-Iter'+str(iter)+'-Num'+str(nmodel)
    mdb.Job(name=jobName, model=modelName, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
        queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, 
        userSubroutine='', scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
        numDomains=numDomainsCPUs, activateLoadBalancing=False, multiprocessingMode=THREADS, numCpus=numDomainsCPUs)     
    mdb.jobs[jobName].submit(consistencyChecking=OFF)
    mdb.jobs[jobName].waitForCompletion()
    print 'The analysis has been completed. On going postprocessing'

    ## Visualization
    myOdb = visualization.openOdb(path = jobName + '.odb')
    myViewport = session.Viewport(name='Viewport-'+modelName, origin =(0.0, 0.0), width=240, height=90)
    myViewport.setValues(displayedObject = myOdb)
    # myViewport.view.setValues(viewOffsetY=offset*nY*ylength)  # Move center of the viewport upward
    myViewport.odbDisplay.display.setValues(plotState=(CONTOURS_ON_UNDEF, CONTOURS_ON_DEF, ))
    myViewport.odbDisplay.commonOptions.setValues(visibleEdges=NONE)
    myViewport.odbDisplay.superimposeOptions.setValues(visibleEdges=NONE)
    myViewport.odbDisplay.commonOptions.setValues(renderStyle = FILLED)
    myViewport.odbDisplay.setPrimaryVariable(variableLabel='U', outputPosition=NODAL, refinement=(COMPONENT, 'U2'))
    myViewport.viewportAnnotationOptions.setValues(triad=OFF, title=OFF, state=OFF, compass=OFF)
    myViewport.viewportAnnotationOptions.setValues(legendFont='-*-verdana-medium-r-normal-*-*-80-*-*-p-*-*-*')
    myViewport.viewportAnnotationOptions.setValues(legendBox=OFF)
    fileName = 'U2-Iter'+str(iter)+'-Num'+str(nmodel)
    session.printToFile(fileName = fileName, format = PNG, canvasObjects = (myViewport, ))
    print 'File ' + fileName + ' has been saved in the working directory. Please check!'

    ## Obtain the last frame of the analysis
    odb = openOdb(path = jobName + '.odb')
    assembly = odb.rootAssembly
    # First step + Last frame of the simulation
    step1 = odb.steps['Step-ThermalExpansion']
    lastFrame = step1.frames[-1]

    UpperCoordinates, LowerCoordinates = [], []
    for i in range(nX):
        nodeName = 'NODE-UPPER-'+str(i)+'-'+str(0)
        nodes = assembly.nodeSets[nodeName]
        displacement = lastFrame.fieldOutputs['U']
        nodesDisplacement = displacement.getSubset(region = nodes)
        # coordinates stores the target coordinates of bottom nodes
        UpperCoordinates.append([UpperXY[i+1][0]+nodesDisplacement.values[0].data[0], 
            UpperXY[i+1][1]+nodesDisplacement.values[0].data[1]])
    for i in range(nXL):
        nodeName = 'NODE-LOWER-'+str(i)+'-'+str(0)
        nodes = assembly.nodeSets[nodeName]
        displacement = lastFrame.fieldOutputs['U']
        nodesDisplacement = displacement.getSubset(region = nodes)
        # coordinates stores the target coordinates of bottom nodes
        LowerCoordinates.append([LowerXY[i+1][0]+nodesDisplacement.values[0].data[0], 
            LowerXY[i+1][1]+nodesDisplacement.values[0].data[1]])

    # Calculate squared error
    def SquaredErrorUpper(xy):
        return (ShapeUpper(xy[0])-xy[1])**2
    def SquaredErrorLower(xy):
        return (ShapeLower(xy[0])-xy[1])**2

    # Defines the target shape
    def ShapeUpper(x):
        a=(1.0*UpperOrigin[1])/(1.0*scaling-UpperOrigin[0])**2
        return a*((x-UpperOrigin[0])**2)+UpperOrigin[1]
    def ShapeLower(x):
        b=(2.0*UpperOrigin[1]-LowerOrigin[1])/(1.0*scaling-LowerOrigin[0])**2
        return b*((x-LowerOrigin[0])**2)+LowerOrigin[1]
    
    RMSE = ((sum(map(SquaredErrorUpper, UpperCoordinates))+sum(map(SquaredErrorLower, LowerCoordinates)))**0.5)/(len(UpperCoordinates)+len(LowerCoordinates))
    del mdb.models[modelName]
    del mdb.jobs[jobName]
    odb.close()
    return RMSE,


creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)
toolbox = base.Toolbox()
## 0: inactive material  1: active material
toolbox.register("attr_bool", random.randint, lower, upper)
# Deformable elements starts building from left to right, bottom to top
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, nX*nY+nXL*nYL)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate", evalOneMax)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
toolbox.register("select", tools.selTournament, tournsize=nsize)

random.seed(64)
pop = toolbox.population(n=sizePopulation)
## Probability of cross and mutation
CXPB, MUTPB = 0.5, 0.2
print("Start of evolution")

### fitness_record contains all the fitness values in the current population at each iteration
fitness_record = []

fitnesses = list(map(toolbox.evaluate, pop, list([0]*sizePopulation), list(range(0,sizePopulation))))
## Store the fitness values of the initial population
fitness_record.append([ind[0] for ind in fitnesses])

for ind, fit in zip(pop, fitnesses):
    ind.fitness.values = fit
print("  Evaluated %i individuals" % len(pop))
fits = [ind.fitness.values[0] for ind in pop]

### Stores the best individual in each iteration
best_ind_all_iterations = []
# Obtain the best individual in the initial population
best_ind_all_iterations.append(list(tools.selBest(pop, 1)[0]))

g = 0
# Begin the evolution
while g < (numIter-1):   # Maximum number of iterations
    # A new generation
    g = g + 1
    print("-- Generation %i --" % g)       
    offspring = toolbox.select(pop, len(pop))
    offspring = list(map(toolbox.clone, offspring))
    # Apply crossover and mutation on the offspring
    for child1, child2 in zip(offspring[::2], offspring[1::2]):
        if random.random() < CXPB:
            # Crossover operation
            toolbox.mate(child1, child2)
            del child1.fitness.values
            del child2.fitness.values
    for mutant in offspring:
        # mutate an individual with probability MUTPB
        if random.random() < MUTPB:
            toolbox.mutate(mutant)
            del mutant.fitness.values
    
    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    # invalid_ind = [ind for ind in offspring]
    fitnesses = map(toolbox.evaluate, invalid_ind, list([g]*len(invalid_ind)), list(range(0,len(invalid_ind))))

    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit
    print("  Evaluated %i individuals" % len(invalid_ind))
    pop[:] = offspring
    ## Stores fitness values for the entire population
    fitnesses_pop = [ind.fitness.values[0] for ind in pop]
    fitness_record.append(fitnesses_pop)
    # Obtain the best individual in the initial population
    best_ind_all_iterations.append(list(tools.selBest(pop, 1)[0]))

print("-- End of (successful) evolution --")
# Obtain the best individual in the current population
best_ind = tools.selBest(pop, 1)[0]
print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values[0]))

with open("fitness_pop.csv", "wb") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(fitness_record)

with open("best_individual.csv", "wb") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['nX','nY','sizePopulation','numIter','expandActive','expandInactive','lower','upper'])
    writer.writerow([nX, nY, sizePopulation, numIter, expandActive, expandInactive, lower, upper])
    writer.writerows(best_ind_all_iterations)
