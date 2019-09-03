#coding:utf-8
from abaqus import *
import testUtils
testUtils.setBackwardCompatibility()
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

### Dimension of a single mesh element
lengthMesh = 1.0
### Number of mesh elements along each dimension of a single deformable element
nx = 15
ny = 1
### Dimensions of a single deformable element
xlength = nx * lengthMesh
ylength = ny * lengthMesh
### Number of deformable elements along each dimension of the global beam
nX = 5
nY = 5
### Create a new model named 'Model-1'
modelName = 'Model-1'
myModel = mdb.Model(name= modelName)
### Create constrained sketch for all deformable elements
for i in range(nX):
    for j in range(nY):
        s = mdb.models[modelName].ConstrainedSketch(name='__profile__', sheetSize=xlength*nX*1.2)
        s.setPrimaryObject(option=STANDALONE)
        s.rectangle(point1=(i*xlength, j*ylength), point2=((i+1)*xlength, (j+1)*ylength))
        # Create part named 'Part-Plate-i-j'
        partName = 'Part-Plate-'+str(i)+'-'+str(j)
        p = mdb.models[modelName].Part(name=partName, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
        p = mdb.models[modelName].parts[partName]
        p.BaseShell(sketch=s)
        s.unsetPrimaryObject()
        # Define vertices sets
        v = p.vertices
        location1 = (i*xlength, j*ylength, 0)
        location2 = (i*xlength, (j+1)*ylength, 0)
        location3 = ((i+1)*xlength, (j+1)*ylength, 0)
        location4 = ((i+1)*xlength, j*ylength, 0)
        verts1 = v.findAt((location1, ))
        verts2 = v.findAt((location2, ))
        verts3 = v.findAt((location3, ))
        verts4 = v.findAt((location4, ))
        p.Set(vertices=verts1, name='Set-1')
        p.Set(vertices=verts2, name='Set-2')
        p.Set(vertices=verts3, name='Set-3')
        p.Set(vertices=verts4, name='Set-4')
        del mdb.models[modelName].sketches['__profile__']


### Create materials named 'Material-Active' and 'Material-Inactive'
activeName = 'Material-Active'
inactiveName = 'Material-Inactive'
mdb.models[modelName].Material(name=activeName)
mdb.models[modelName].materials[activeName].Elastic(table=((210000.0, 0.3), ))
mdb.models[modelName].materials[activeName].Expansion(table=((0.3, ), ))
mdb.models[modelName].Material(name=inactiveName)
mdb.models[modelName].materials[inactiveName].Elastic(table=((210000.0, 0.3), ))
mdb.models[modelName].materials[inactiveName].Expansion(table=((0.1, ), ))
### Create HomogeneousSolidSection objects and assign material property to them 
activeSection = 'Section-Active'
inactiveSection = 'Section-Inactive'
mdb.models[modelName].HomogeneousSolidSection(name=activeSection, material=activeName, thickness=None)
mdb.models[modelName].HomogeneousSolidSection(name=inactiveSection, material=inactiveName, thickness=None)


### Adding sections to parts accordingly
### Resulting section pattern is like a checkerboard
for i in range(nX):
    for j in range(nY):
        # Create partName named 'Part-Plate-i-j'
        partName = 'Part-Plate-'+str(i)+'-'+str(j)
        p = mdb.models[modelName].parts[partName]
        f = p.faces
        faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
        # Create a region
        region = regionToolset.Region(faces=faces)
        if (i+j)%2==0:
            p.SectionAssignment(region=region, sectionName=inactiveSection, 
                offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
        else:
            p.SectionAssignment(region=region, sectionName=activeSection, 
                offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)


### Create assembly consisting of all plate instances
a = mdb.models[modelName].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
for i in range(nX):
    for j in range(nY):
        # Create partName named 'Part-Plate-i-j'
        partName = 'Part-Plate-'+str(i)+'-'+str(j)
        instanceName = partName+'-1'
        p = mdb.models[modelName].parts[partName]
        # Create indenpendent instances by setting 'dependent' to OFF
        a.Instance(name=instanceName, part=p, dependent=OFF)


### Generate mesh for the assembly
# Create tuple of part instances
partInstances = []
for i in range(nX):
    for j in range(nY):
        instanceName = 'Part-Plate-'+str(i)+'-'+str(j)+'-1'
        partInstances.append(a.instances[instanceName])
partInstances=tuple(partInstances)
# Seed part instances
a.seedPartInstance(regions=partInstances, size=lengthMesh, deviationFactor=0.1, minSizeFactor=0.1)
# Generate mesh for partInstances based on previous seeding
a.generateMesh(regions=partInstances)


### Create a new step named 'Step-ThermalExpansion' after initial step
mdb.models[modelName].StaticStep(name='Step-ThermalExpansion', previous='Initial')


### Create interaction edges and define interaction relations 
for i in range(nX):
    for j in range(nY):
        # Create instanceName named 'Part-Plate-i-j-1'
        instanceName = 'Part-Plate-'+str(i)+'-'+str(j)+'-1'
        s1 = a.instances[instanceName].edges
        edgeCenter1 = ((i+0.5)*xlength, (j+1.0)*ylength, 0)
        edgeCenter2 = ((i+0.0)*xlength, (j+0.5)*ylength, 0)
        edgeCenter3 = ((i+1.0)*xlength, (j+0.5)*ylength, 0)
        edgeCenter4 = ((i+0.5)*xlength, (j+0.0)*ylength, 0)
        edges1 = s1.findAt((edgeCenter1, ))
        edges2 = s1.findAt((edgeCenter2, ))
        edges3 = s1.findAt((edgeCenter3, ))
        edges4 = s1.findAt((edgeCenter4, ))
        edges1 = (edges1, ); edges2 = (edges2, ); edges3 = (edges3, ); edges4 = (edges4, )
        # Create surfaces for all 4 edges of current element
        # a is rootAssembly
        a.Surface(side1Edges=edges1, name='Surf-Plate-'+str(i)+'-'+str(j)+'-1')
        a.Surface(side1Edges=edges2, name='Surf-Plate-'+str(i)+'-'+str(j)+'-2')
        a.Surface(side1Edges=edges3, name='Surf-Plate-'+str(i)+'-'+str(j)+'-3')
        a.Surface(side1Edges=edges4, name='Surf-Plate-'+str(i)+'-'+str(j)+'-4')
# Define interaction relations between lower and upper surfaces
for i in range(nX):
    for j in range(nY-1):
        surfName1 = 'Surf-Plate-'+str(i)+'-'+str(j)+'-1'
        surfName2 = 'Surf-Plate-'+str(i)+'-'+str(j+1)+'-4'
        region1=a.surfaces[surfName1]
        region2=a.surfaces[surfName2]
        contactName = 'Constraint-'+str(i)+str(j)+str(i)+str(j+1)+'-14'
        mdb.models[modelName].Tie(name=contactName, master=region1, slave=region2, 
            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
# Define interaction relations between left and right surfaces
for i in range(nX-1):
    for j in range(nY):
        surfName1 = 'Surf-Plate-'+str(i)+'-'+str(j)+'-3'
        surfName2 = 'Surf-Plate-'+str(i+1)+'-'+str(j)+'-2'
        region1=a.surfaces[surfName1]
        region2=a.surfaces[surfName2]
        contactName = 'Constraint-'+str(i)+str(j)+str(i+1)+str(j)+'-32'
        mdb.models[modelName].Tie(name=contactName, master=region1, slave=region2, 
            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)


### Create boundary conditions for the assembly
a = mdb.models[modelName].rootAssembly
for j in range(nY):
    instanceName = 'Part-Plate-0-'+str(j)+'-1'
    e1 = a.instances[instanceName].edges
    edgeCenter = (0.0, (j+0.5)*ylength, 0.0)
    edges1 = e1.findAt((edgeCenter, ))
    edges1 = (edges1, )
    setName = 'Set-Plate-0-'+str(j)+'-Fix'
    a.Set(edges=edges1, name=setName)
    # Set boundary conditions
    region = a.sets[setName]
    bcName = 'BC-Plate-0-'+str(j)+'-Fix'
    mdb.models[modelName].DisplacementBC(name=bcName, 
        createStepName='Initial', region=region, u1=SET, u2=UNSET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
instanceName = 'Part-Plate-0-0-1'
v1 = a.instances[instanceName].vertices
vertexCenter = (0.0, 0.0, 0.0)
v = v1.findAt((vertexCenter, ))
v = (v, )
setName = 'Set-Plate-Origin-Fix'
a.Set(vertices=v, name=setName)
region = a.sets[setName]
bcName = 'BC-Origin-Fix'
# Origin fix in x-direction, y-direction, xy-rotation
mdb.models[modelName].DisplacementBC(name=bcName, 
        createStepName='Initial', region=region, u1=SET, u2=SET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)


### Define nodeSets that records the displacement information of the entire assembly
# 1st: bottom left vertices of all deformable elements
a = mdb.models[modelName].rootAssembly
for i in range(nX):
    for j in range(nY):
        nodeName = 'Node-Plate-'+str(i)+'-'+str(j)
        instanceName = 'Part-Plate-'+str(i)+'-'+str(j)+'-1'
        v1 = a.instances[instanceName].vertices
        vertexCenter = (i*xlength, j*ylength, 0.0)
        v = v1.findAt((vertexCenter, ))
        v = (v, )
        a.Set(vertices=v, name=nodeName)
# 2nd: top left vertices of elements in the top row
for i in range(nX):
        nodeName = 'Node-Plate-'+str(i)+'-'+str(nY)
        instanceName = 'Part-Plate-'+str(i)+'-'+str(nY-1)+'-1'
        v1 = a.instances[instanceName].vertices
        vertexCenter = (i*xlength, nY*ylength, 0.0)
        v = v1.findAt((vertexCenter, ))
        v = (v, )
        a.Set(vertices=v, name=nodeName)
# 3rd: bottom right vertices of elements in the rightmost column
for j in range(nY):
        nodeName = 'Node-Plate-'+str(nX)+'-'+str(j)
        instanceName = 'Part-Plate-'+str(nX-1)+'-'+str(j)+'-1'
        v1 = a.instances[instanceName].vertices
        vertexCenter = (nX*xlength, j*ylength, 0.0)
        v = v1.findAt((vertexCenter, ))
        v = (v, )
        a.Set(vertices=v, name=nodeName)
# 4th: top right vertex of the entire assembly
nodeName = 'Node-Plate-'+str(nX)+'-'+str(nY)
instanceName = 'Part-Plate-'+str(nX-1)+'-'+str(nY-1)+'-1'
v1 = a.instances[instanceName].vertices
vertexCenter = (nX*xlength, nY*ylength, 0.0)
v = v1.findAt((vertexCenter, ))
v = (v, )
a.Set(vertices=v, name=nodeName)


### Define ambient temperature field
a = mdb.models[modelName].rootAssembly
allFaceSet = []
for i in range(nX):
    for j in range(nY):
        instanceName = 'Part-Plate-'+str(i)+'-'+str(j)+'-1'
        f1 = a.instances[instanceName].faces
        faces = f1.getSequenceFromMask(mask=('[#1 ]', ), )
        allFaceSet.append(faces)
allFaceSet =tuple(allFaceSet)
setName = 'Set-allFace'
a.Set(faces=allFaceSet, name=setName)
region = a.sets[setName]
# Set temperature for initial and final state
initialT = 20.0
finalT = 120.0
mdb.models[modelName].Temperature(name='Predefined Field-Temperature', 
    createStepName='Initial', region=region, distributionType=UNIFORM, 
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(initialT, ))
mdb.models[modelName].predefinedFields['Predefined Field-Temperature'].setValuesInStep(
    stepName='Step-ThermalExpansion', magnitudes=(finalT, ))


### Create job and submit
jobName = 'Job-MultiPlateContact'
mdb.Job(name=jobName, model=modelName, description='', 
    type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
    memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', resultsFormat=ODB)
mdb.jobs[jobName].submit(consistencyChecking=OFF)
mdb.jobs[jobName].waitForCompletion()
print 'The analysis has been completed. On going postprocessing'


### Visualization
myOdb = visualization.openOdb(path = jobName + '.odb')
myViewport = session.Viewport(name='Viewport-'+modelName, origin =(0.5*nX*xlength, 0.5*nY*ylength), 
    width=240, height=120)
myViewport.setValues(displayedObject = myOdb)
myViewport.odbDisplay.display.setValues(plotState = CONTOURS_ON_DEF)
myViewport.odbDisplay.commonOptions.setValues(renderStyle = FILLED)
session.printToFile(fileName = 'Mises', format = PNG, canvasObjects = (myViewport, ))
print 'File Mises.png has been saved in the working directory. Please check!'



